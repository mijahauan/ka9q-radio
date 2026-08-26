
// Version that just forks and execs the commands 'avahi-browse' or 'avahi-publish'
// Using the Avahi API is a nightmare; it's a twisted mess of callbacks between "callbacks", "services", "resolvers", etc.
// March 2024, Phil Karn KA9Q
// Jan 2025 - option to publish through static additions to avahi daemon files /etc/avahi/hosts, /etc/avahi/services
//            note - does not remove such entries yet!
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <arpa/inet.h>
#include <string.h>
#include <errno.h>
#include <sys/file.h>
#include <stdbool.h>
#include <ctype.h>
#include "misc.h"
#include "avahi.h"

#define SERVICES "/etc/avahi/services"
#define HOSTS "/etc/avahi/hosts"

bool Static_avahi;

#ifdef __APPLE__
// macOS has no avahi; register with the system mDNSResponder through the
// native dns_sd API instead of exec'ing the avahi CLI tools.
// Registrations must stay active for the life of the process (like the
// forked avahi-publish children on Linux), so each avahi_start() hands its
// DNSServiceRefs to a detached thread running their event loop and nothing
// is ever deallocated. The Static_avahi (/etc/avahi) mode doesn't apply here.
#include <dns_sd.h>
#include <pthread.h>
#include <sys/select.h>
#include <netinet/in.h>

struct dnssd_session {
  DNSServiceRef refs[2]; // service registration + connection holding the address record
  int nrefs;
};

static void *dnssd_event_loop(void *arg){
  struct dnssd_session * const sp = arg;
  pthread_setname("dns-sd");
  while(true){
    fd_set rfds;
    FD_ZERO(&rfds);
    int maxfd = -1;
    for(int i = 0; i < sp->nrefs; i++){
      int const fd = DNSServiceRefSockFD(sp->refs[i]);
      FD_SET(fd,&rfds);
      if(fd > maxfd)
	maxfd = fd;
    }
    if(select(maxfd+1,&rfds,NULL,NULL,NULL) <= 0)
      return NULL;
    for(int i = 0; i < sp->nrefs; i++){
      if(FD_ISSET(DNSServiceRefSockFD(sp->refs[i]),&rfds)
	 && DNSServiceProcessResult(sp->refs[i]) != kDNSServiceErr_NoError)
	return NULL;
    }
  }
}

static void dnssd_register_callback(DNSServiceRef sdref,DNSServiceFlags flags,DNSServiceErrorType error,
				    char const *name,char const *regtype,char const *domain,void *context){
  (void)sdref; (void)flags; (void)context;
  if(error != kDNSServiceErr_NoError)
    fprintf(stderr,"dns-sd registration of %s.%s%s failed: error %d\n",name,regtype,domain,(int)error);
}

static void dnssd_record_callback(DNSServiceRef sdref,DNSRecordRef rref,DNSServiceFlags flags,
				  DNSServiceErrorType error,void *context){
  (void)sdref; (void)rref; (void)flags;
  if(error != kDNSServiceErr_NoError)
    fprintf(stderr,"dns-sd address record for %s failed: error %d\n",(char *)context,(int)error);
}

// Append a string to a TXT record as a length-prefixed entry
static int txt_append(uint8_t *txt,int len,int size,char const *s){
  int const slen = strlen(s);
  if(slen == 0 || slen > 255 || len + 1 + slen > size)
    return len;
  txt[len] = slen;
  memcpy(txt + len + 1,s,slen);
  return len + 1 + slen;
}

int avahi_start(char const *service_name,char const *service_type,int const service_port,char const *dns_name,uint32_t address,char const *description){
  // Same TXT entries the Linux path passes to avahi-publish-service
  uint8_t txt[512];
  int txtlen = 0;
  if(description != NULL)
    txtlen = txt_append(txt,txtlen,sizeof(txt),description);
  char tmp[256];
  snprintf(tmp,sizeof(tmp),"pid=%d",getpid());
  txtlen = txt_append(txt,txtlen,sizeof(txt),tmp);
  char hostname[sysconf(_SC_HOST_NAME_MAX)];
  gethostname(hostname,sizeof(hostname));
  snprintf(tmp,sizeof(tmp),"source=%s",hostname);
  txtlen = txt_append(txt,txtlen,sizeof(txt),tmp);

  struct dnssd_session * const sp = calloc(1,sizeof(*sp));
  DNSServiceErrorType err = DNSServiceRegister(&sp->refs[sp->nrefs],
					       0, // default flags: rename on conflict
					       kDNSServiceInterfaceIndexAny,
					       service_name,
					       service_type,
					       NULL, // domain: local
					       dns_name, // SRV target, like avahi-publish --host
					       htons(service_port),
					       txtlen,txt,
					       dnssd_register_callback,NULL);
  if(err != kDNSServiceErr_NoError)
    fprintf(stderr,"DNSServiceRegister(%s,%s) failed: error %d\n",service_name,service_type,(int)err);
  else
    sp->nrefs++;

  if(address != 0){
    // Publish an address record mapping dns_name to the (multicast) address,
    // like avahi-publish-address on Linux
    DNSServiceRef conn = NULL;
    err = DNSServiceCreateConnection(&conn);
    if(err != kDNSServiceErr_NoError){
      fprintf(stderr,"DNSServiceCreateConnection failed: error %d\n",(int)err);
    } else {
      uint32_t const addr_be = htonl(address);
      DNSRecordRef rref = NULL;
      err = DNSServiceRegisterRecord(conn,&rref,
				     kDNSServiceFlagsUnique,
				     kDNSServiceInterfaceIndexAny,
				     dns_name,
				     kDNSServiceType_A,
				     kDNSServiceClass_IN,
				     sizeof(addr_be),&addr_be,
				     0, // default TTL
				     dnssd_record_callback,strdup(dns_name));
      if(err != kDNSServiceErr_NoError){
	fprintf(stderr,"DNSServiceRegisterRecord(%s) failed: error %d\n",dns_name,(int)err);
	DNSServiceRefDeallocate(conn);
      } else
	sp->refs[sp->nrefs++] = conn;
    }
  }
  if(sp->nrefs == 0){
    free(sp);
    return -1;
  }
  pthread_t thread;
  pthread_create(&thread,NULL,dnssd_event_loop,sp);
  pthread_detach(thread);
  return 0;
}

#else // !__APPLE__: fork/exec the avahi CLI tools, or write static /etc/avahi files

int avahi_publish_address(char const *dns_name,char const *ip_address_string);
int avahi_publish_service(char const *service_name, char const *service_type, char const *dns_name,int service_port, char const *description, int pid);
static void sanitize_filename(char *name);

int avahi_start(char const *service_name,char const *service_type,int const service_port,char const *dns_name,uint32_t address,char const *description){
  int pid = getpid(); // Advertise the parent's pid, not the child's
  if(Static_avahi){
    // Just put in /etc/avahi/services
    avahi_publish_service(service_name, service_type, dns_name,service_port,description,pid);
  } else { // static = no
    if(fork() == 0){
#if 0
      fprintf(stderr,"avahi-publish-service child pid %d\n",getpid());
#endif
      // run "avahi-publish-service --no-fail --host=dns_name service_name service_type service_port description pid hostname" in subprocess
      // No need to free the asprintf-allocated strings, we're calling exec anyway
      // Turn off warnings about unused return values from asprintf generated by some compilers
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
      char *port_string = NULL;
      asprintf(&port_string, "%d", service_port);

      char *host_string = NULL;
      asprintf(&host_string, "--host=%s",dns_name);

      char *pid_string = NULL;
      asprintf(&pid_string, "pid=%d", pid);

      char hostname[sysconf(_SC_HOST_NAME_MAX)];
      gethostname(hostname,sizeof(hostname));
      char *source_string = NULL;
      asprintf(&source_string, "source=%s", hostname);
#pragma GCC diagnostic pop

#if 1
      fprintf(stderr,"%s %s %s %s %s %s %s %s %s %s\n",
	      "avahi-publish-service", "avahi-publish-service", "--no-fail", host_string, service_name, service_type, port_string, description, pid_string, hostname);
#endif
      execlp("avahi-publish-service", "avahi-publish-service", "--no-fail", host_string, service_name, service_type, port_string, description, pid_string, source_string, NULL);
      perror("exec avahi publish service failed");
      return -1;
    }
  }
  if(address != 0){
    char *ip_address_string = NULL;    // No need to free, we're calling exec
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
    asprintf(&ip_address_string,"%d.%d.%d.%d",(address >> 24) & 0xff, (address >> 16) & 0xff, (address >> 8) & 0xff, address & 0xff);
#pragma GCC diagnostic pop
    if(Static_avahi){
      // Just put in /etc/avahi/hosts
      avahi_publish_address(dns_name,ip_address_string);
    } else {
      if(fork() == 0){
	// run "avahi-publish-address dns_name address", only if an address is specified

#if 0
	fprintf(stderr,"avahi start: ip address string = %s\n",ip_address_string);
	fprintf(stderr,"avahi-publish-address child pid %d\n",getpid());
	fprintf(stderr,"%s %s %s %s\n",
		"avahi-publish-address", "avahi-publish-address",dns_name,ip_address_string);
#endif
	execlp("avahi-publish-address", "avahi-publish-address",dns_name,ip_address_string, NULL);
	perror("exec avahi publish address");
	return -1;
      }
    }
  }
  return 0;
}

// Publish service in static avahi directory
int avahi_publish_service(char const *service_name, char const *service_type, char const *dns_name,int service_port, char const *description, int pid){

  char *sanitized_name = strdup(service_name);
  sanitize_filename(sanitized_name);
  char *sanitized_type = strdup(service_type);
  sanitize_filename(sanitized_type);

  char tmp_file[1024];
  snprintf(tmp_file,sizeof(tmp_file),"%s/%s-%s.service.tmp%d",SERVICES,sanitized_name,sanitized_type,pid);
  FILE *fp = fopen(tmp_file,"w"); // Will overwrite if exists
  if(fp == NULL){
    fprintf(stderr,"Can't create %s: %s\n",tmp_file,strerror(errno));
    FREE(sanitized_type);
    FREE(sanitized_name);
    return -1;
  }
  char hostname[sysconf(_SC_HOST_NAME_MAX)];
  gethostname(hostname,sizeof(hostname));

  fputs("<service-group>\n",fp);
  fprintf(fp,"<name>%s</name>\n",service_name);
  fputs("<service protocol=\"ipv4\">\n",fp);
  fprintf(fp,"<host-name>%s</host-name>\n",dns_name);
  fprintf(fp,"<type>%s</type>\n",service_type);
  fprintf(fp,"<port>%d</port>\n",service_port);
  fprintf(fp,"<txt-record>pid=%d</txt-record>\n",pid);
  fprintf(fp,"<txt-record>%s</txt-record>\n",description);
  fprintf(fp,"<txt-record>source=%s</txt-record>\n",hostname);
  fputs("</service>\n",fp);
  fputs("</service-group>\n",fp);
  fclose(fp);
  char service_file[1024];
  snprintf(service_file,sizeof(service_file),"%s/%s-%s.service",SERVICES,sanitized_name,sanitized_type);

  int r = rename(tmp_file,service_file);
  if(r != 0)
    fprintf(stderr,"Can't rename %s to %s: %s\n",
	    tmp_file,service_file,strerror(errno));
  FREE(sanitized_type);
  FREE(sanitized_name);

  return r;
}

// Publish address, host pair in static avahi file

int avahi_publish_address(char const *name,char const *address){
  FILE *fp = fopen(HOSTS,"a+");
  if(fp == NULL){
    fprintf(stderr,"Can't open %s: %s\n",HOSTS,strerror(errno));
    return -1;
  }
  // the lock has to be exclusive even when reading, otherwise there could be double entries
  if(flock(fileno(fp),LOCK_EX) != 0){
    fprintf(stderr,"Can't lock %s: %s\n",HOSTS,strerror(errno));
    fclose(fp);
    return -1;
  }
  rewind(fp); // Start reading at front
  char *buffer = NULL;
  size_t linecap = 0;
  ssize_t linelen;
  while((linelen = getline(&buffer,&linecap,fp)) > 0){
    char *cp;
    if((cp = strchr(buffer,'\n')) != NULL)
      *cp = '\0';
    if(buffer[0] == '#')
      continue;
    char *line = buffer;
    char const *faddr = strsep(&line," \t");
    if(faddr == NULL || line == NULL)
      continue;

    // Skip multiple white space
    while(*line == ' ' || *line == '\t')
      line++;

    char const *fname = line;
    if(fname == NULL || *fname == '\0')
      continue;

    // These are string comparisons so there might be a false mismatch
    // if they're written differently. redo as compares of canonical forms
    // Also, what are the rules about the same names with different addresses
    // and vice versa?
    if(strcmp(name,fname) == 0 && strcmp(address,faddr) == 0){
      // Yes, entry is already present, don't do anything
      flock(fileno(fp),LOCK_UN);
      fclose(fp);
      free(buffer);
      return 0;
    }
  }
  free(buffer);
  // we've read the entire file without a match, so append our record
  // The "a+" open mode means we'll append
  if(fprintf(fp,"%s %s\n",address,name) < 0)
    fprintf(stderr,"Can't append to %s: %s\n",HOSTS,strerror(errno));
  flock(fileno(fp),LOCK_UN); // Upgrade to exclusive lock for writing
  fclose(fp);

  return 0;
}
// Remove problematic characters from an Avahi service file name
static void sanitize_filename(char *name) {
  char *src = name;
  char *dst = name;

  // Trim leading whitespace
  while (isspace((unsigned char)*src))
    src++;

  while (*src) {
    unsigned char c = (unsigned char)*src;
    if (isalnum(c)) {
      *dst++ = c;
    } else if (c == ' ' || c == '_' || c == '-' || c == '@') {
      *dst++ = '-';  // Normalize all to hyphen
    } else if (c == '.') {
      *dst++ = '.';  // Keep dot for extension
    } else {
      // Drop anything else (e.g., control chars, punctuation)
    }
    src++;
  }
  // Trim trailing hyphens or whitespace
  while (dst > name && (dst[-1] == '-' || isspace((unsigned char)dst[-1])))
    dst--;
  *dst = '\0';
}
#endif // !__APPLE__
