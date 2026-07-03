#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "misc.h"
#include "avahi.h"


#if 0 // devel testing
int main(){
  char const *service_name = "_ka9q-ctl._udp";

  struct service_tab table[1000];
  int line_count = avahi_browse(table,1000,service_name);


  for(int i=0;i<line_count;i++){
    struct service_tab *tp = &table[i];
    printf("line %s interface %s protocol %s name %s type %s domain %s dns %s address %s port %s txt %s\n",
	   tp->line_type,
	   tp->interface,
	   tp->protocol,
	   tp->name,
	   tp->type,
	   tp->domain,
	   tp->dns_name,
	   tp->address,
	   tp->port,
	   tp->txt);
  }
  exit(0);
}
#endif

#ifdef __APPLE__
// The system responder reports each service on every interface, including
// virtual ones (feth, awdl, utun, bridge) that can't be usefully joined.
// Rank physical interfaces (en*) first so deduplication keeps one of those.
static int iface_rank(char const *ifname){
  return (ifname != NULL && ifname[0] == 'e' && ifname[1] == 'n') ? 0 : 1;
}
#endif

// NULL entries sort to end so we can drop them after sorting and removing dupes
static int table_compare(void const *a,void const *b){
  struct service_tab const *t1 = a;
  struct service_tab const *t2 = b;
  if(t2 == NULL || t2->name == NULL)
    return -1;
  else if(t1 == NULL || t1->name == NULL)
    return +1;
  int const r = strcmp(t1->name,t2->name);
#ifdef __APPLE__
  if(r == 0)
    return iface_rank(t1->interface) - iface_rank(t2->interface);
#endif
  return r;
}

// Sort by instance name, remove duplicates (same service seen on multiple
// interfaces), return the number of surviving entries
static int sort_and_dedupe(struct service_tab *table,int line_count){
  qsort(table,line_count,sizeof(table[0]),table_compare);

  int dupes = 0;
  for(int i = 0; i < line_count; i++){
    if(table[i].buffer == NULL || table[i].name == NULL) // already zeroed out?
      continue;
    for(int j = i+1; j < line_count; j++){
      if(table[j].buffer != NULL && table[j].name != NULL){
	if(strcmp(table[i].name,table[j].name) != 0)
	  break; // Different; stop looking
	// zero it out, count it and look for more
	FREE(table[j].buffer);
	memset(&table[j], 0, sizeof table[j]);
	dupes++;
      }
    }
  }
  // re-sort, pushing zeroed entries to end
  qsort(table,line_count,sizeof(table[0]),table_compare);
  return line_count - dupes;
}
#ifdef __APPLE__
// macOS has no avahi-browse; query the system mDNSResponder directly
// through the native dns_sd API. All operations share one connection so a
// single-threaded DNSServiceProcessResult() loop dispatches every callback:
// browse results start resolves, resolves start address queries, and
// completed address queries fill the caller's table with the same fields
// the avahi-browse output parser produces on Linux.
#include <dns_sd.h>
#include <net/if.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <time.h>
#include <unistd.h>

struct browse_state {
  struct service_tab *table;
  int tabsize;
  int count;
  DNSServiceRef connection;
  int pending;               // outstanding resolve/getaddrinfo operations
  struct resolve_state *ops; // all per-service operations, freed at the end
};

struct resolve_state {
  struct resolve_state *next;
  struct browse_state *bs;
  DNSServiceRef ref;     // subordinate ref sharing bs->connection
  DNSServiceRef gai_ref; // ditto, for the address query
  uint32_t interface_index;
  bool resolved;
  bool recorded;
  char *name;
  char *type;
  char *domain;
  char *dns_name;
  uint16_t port; // host byte order
};

static void gai_callback(DNSServiceRef sdref,DNSServiceFlags flags,uint32_t ifindex,
			 DNSServiceErrorType error,char const *hostname,
			 struct sockaddr const *sa,uint32_t ttl,void *context){
  (void)sdref; (void)flags; (void)hostname; (void)ttl;
  struct resolve_state * const rs = context;
  struct browse_state * const bs = rs->bs;
  if(error == kDNSServiceErr_NoError && !rs->recorded && sa != NULL && sa->sa_family == AF_INET
     && bs->count < bs->tabsize){
    struct sockaddr_in const *sin = (struct sockaddr_in const *)sa;
    char address[INET_ADDRSTRLEN];
    inet_ntop(AF_INET,&sin->sin_addr,address,sizeof(address));
    char port[16];
    snprintf(port,sizeof(port),"%u",rs->port);
    char ifname[IF_NAMESIZE] = "";
    uint32_t const ifx = ifindex != 0 ? ifindex : rs->interface_index;
    if(ifx == 0 || if_indextoname(ifx,ifname) == NULL)
      ifname[0] = '\0';

    // One buffer holds all the strings, like the getline buffer on Linux,
    // so avahi_free_service_table() works unchanged
    struct service_tab * const tp = &bs->table[bs->count];
    size_t const need = sizeof("=") + sizeof("IPv4") + strlen(ifname) + strlen(rs->name)
      + strlen(rs->type) + strlen(rs->domain) + strlen(rs->dns_name)
      + strlen(address) + strlen(port) + 8; // separating NULs (txt is empty)
    char *wp = malloc(need);
    tp->buffer = wp;
    tp->line_type = wp; wp = stpcpy(wp,"=") + 1;
    tp->interface = wp; wp = stpcpy(wp,ifname) + 1;
    tp->protocol = wp;  wp = stpcpy(wp,"IPv4") + 1;
    tp->name = wp;      wp = stpcpy(wp,rs->name) + 1;
    tp->type = wp;      wp = stpcpy(wp,rs->type) + 1;
    tp->domain = wp;    wp = stpcpy(wp,rs->domain) + 1;
    tp->dns_name = wp;  wp = stpcpy(wp,rs->dns_name) + 1;
    tp->address = wp;   wp = stpcpy(wp,address) + 1;
    tp->port = wp;      wp = stpcpy(wp,port) + 1;
    tp->txt = wp;       *wp = '\0';
    bs->count++;
    rs->recorded = true;
    bs->pending--;
  } else if(error != kDNSServiceErr_NoError && !rs->recorded){
    rs->recorded = true; // no address is coming; stop waiting for this one
    bs->pending--;
  }
}

static void resolve_callback(DNSServiceRef sdref,DNSServiceFlags flags,uint32_t ifindex,
			     DNSServiceErrorType error,char const *fullname,char const *hosttarget,
			     uint16_t port_be,uint16_t txtlen,unsigned char const *txt,void *context){
  (void)sdref; (void)flags; (void)fullname; (void)txtlen; (void)txt;
  struct resolve_state * const rs = context;
  struct browse_state * const bs = rs->bs;
  if(rs->resolved)
    return; // duplicate callback
  if(error != kDNSServiceErr_NoError){
    rs->resolved = rs->recorded = true;
    bs->pending--;
    return;
  }
  rs->resolved = true;
  rs->dns_name = strdup(hosttarget);
  rs->port = ntohs(port_be);
  rs->gai_ref = bs->connection;
  if(DNSServiceGetAddrInfo(&rs->gai_ref,kDNSServiceFlagsShareConnection,ifindex,
			   kDNSServiceProtocol_IPv4,hosttarget,gai_callback,rs) != kDNSServiceErr_NoError){
    rs->gai_ref = NULL;
    rs->recorded = true;
    bs->pending--;
  }
}

static void browse_callback(DNSServiceRef sdref,DNSServiceFlags flags,uint32_t ifindex,
			    DNSServiceErrorType error,char const *name,char const *regtype,
			    char const *domain,void *context){
  (void)sdref;
  struct browse_state * const bs = context;
  if(error != kDNSServiceErr_NoError || !(flags & kDNSServiceFlagsAdd))
    return;
  struct resolve_state * const rs = calloc(1,sizeof(*rs));
  rs->bs = bs;
  rs->interface_index = ifindex;
  rs->name = strdup(name);
  rs->type = strdup(regtype);
  rs->domain = strdup(domain);
  rs->ref = bs->connection;
  if(DNSServiceResolve(&rs->ref,kDNSServiceFlagsShareConnection,ifindex,
		       name,regtype,domain,resolve_callback,rs) == kDNSServiceErr_NoError){
    bs->pending++;
    rs->next = bs->ops;
    bs->ops = rs;
  } else {
    free(rs->name);
    free(rs->type);
    free(rs->domain);
    free(rs);
  }
}

int avahi_browse(struct service_tab *table,int tabsize,char const *service_name){
  if(service_name == NULL || strlen(service_name) == 0 || table == NULL || tabsize == 0)
    return 0;

  memset(table,0,sizeof(*table) * tabsize); // Clear it out in case we free it later

  struct browse_state bs;
  memset(&bs,0,sizeof(bs));
  bs.table = table;
  bs.tabsize = tabsize;
  if(DNSServiceCreateConnection(&bs.connection) != kDNSServiceErr_NoError)
    return 0;
  DNSServiceRef browse = bs.connection;
  if(DNSServiceBrowse(&browse,kDNSServiceFlagsShareConnection,kDNSServiceInterfaceIndexAny,
		      service_name,NULL,browse_callback,&bs) != kDNSServiceErr_NoError){
    DNSServiceRefDeallocate(bs.connection);
    return 0;
  }
  // Collect responses. mDNS answers usually arrive within a few hundred ms;
  // leave once everything found has resolved, with a hard cap comparable to
  // an avahi-browse -t run.
  int const fd = DNSServiceRefSockFD(bs.connection);
  struct timespec start;
  clock_gettime(CLOCK_MONOTONIC,&start);
  while(bs.count < tabsize){
    struct timespec now;
    clock_gettime(CLOCK_MONOTONIC,&now);
    long const elapsed_ms = (now.tv_sec - start.tv_sec) * 1000 + (now.tv_nsec - start.tv_nsec) / 1000000;
    if(elapsed_ms >= 2500)
      break;
    if(elapsed_ms >= 1000 && bs.pending == 0 && bs.count > 0)
      break;
    fd_set rfds;
    FD_ZERO(&rfds);
    FD_SET(fd,&rfds);
    struct timeval timeout = { .tv_sec = 0, .tv_usec = 100000 };
    int const n = select(fd+1,&rfds,NULL,NULL,&timeout);
    if(n < 0)
      break;
    if(n > 0 && DNSServiceProcessResult(bs.connection) != kDNSServiceErr_NoError)
      break;
  }
  DNSServiceRefDeallocate(bs.connection); // also cancels all operations sharing it
  for(struct resolve_state *rs = bs.ops; rs != NULL;){
    struct resolve_state * const next = rs->next;
    free(rs->name);
    free(rs->type);
    free(rs->domain);
    free(rs->dns_name);
    free(rs);
    rs = next;
  }
  return sort_and_dedupe(table,bs.count);
}

#else // !__APPLE__: parse the output of the avahi-browse CLI tool

// De-escape decimal sequences of the form \032 to (space)
// Do it in place since the string can only get shorter
static int deescape(char *s){
  int len = 0;
  char *rp = s;
  char *wp = s;
  bool add_null = false;

  while(*rp != '\0'){
    if(*rp == '\\'){
      rp++;
      *wp = (*rp++ - '0') * 100;
      *wp += (*rp++ - '0') * 10;
      *wp++ += (*rp++ - '0');
      add_null = true;
    } else {
      *wp++ = *rp++; // Copy original character, which may have shifted
    }
    len++;
  }
  if(add_null)
    *wp = '\0'; // String has shortened, make sure it's null terminated
  return len;
}



// Invoke 'avahi-browse' command, filter and sort output
int avahi_browse(struct service_tab *table,int tabsize,char const *service_name){
  if(service_name == NULL || strlen(service_name) == 0 || table == NULL || tabsize == 0)
    return 0;

  memset(table,0,sizeof(*table) * tabsize); // Clear it out in case we free it later

  FILE *fp;
  {
    char cmd[1024];
    snprintf(cmd,sizeof(cmd),"avahi-browse -ptr %s",service_name);
    fp = popen(cmd,"r");
  }
  if(fp == NULL)
    exit(1);

  int line_count;
  char *line = NULL;
  size_t linesize = 0;

  for(line_count = 0; line_count < tabsize;){
    struct service_tab *tp = &table[line_count];

    // Allocates or reallocates as necessary
    if(getline(&line,&linesize,fp) <= 0){
      FREE(line);
      break;
    }
    tp->buffer = line; // so it can be freed later
    tp->line_type = strsep(&line,";");
    tp->interface = strsep(&line,";");
    tp->protocol =  strsep(&line,";");
    tp->name =  strsep(&line,";");
    tp->type =  strsep(&line,";");
    tp->domain = strsep(&line,";");
    tp->dns_name = strsep(&line,";");
    tp->address = strsep(&line,";");
    tp->port =  strsep(&line,";");
    tp->txt =   strsep(&line,";");

    if(strcmp(tp->line_type,"=") == 0 // Only full lines with resolved addresses
       && strcmp(tp->protocol,"IPv4") == 0 // Only IPv4 for now
       && strcmp(tp->type,service_name) == 0){ // Should always match, but to be sure
      // Keep this line and the pointers into it
      deescape(tp->name);
      tp++;
      line_count++;
      line = NULL; // Force a new allocation by getline() on next iteration
      linesize = 0;
    }
    // Otherwise reuse line (with possible expansion)
  }
  FREE(line);
  pclose(fp); // What to do with return code?
  return sort_and_dedupe(table,line_count);
}
#endif // !__APPLE__

// Free and wipe pointers inside service table; caller must free table array itself
void avahi_free_service_table(struct service_tab *table,int tabsize){
  for(int i=0; i < tabsize; i++)
    free(table[i].buffer);

  memset(table,0,sizeof(*table) * tabsize); // Clear it all out
}
