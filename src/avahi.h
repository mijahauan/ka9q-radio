#ifndef _AVAHI_H
// The publish/browse implementations shell out to the avahi CLI tools on
// Linux and use the native dns_sd API on macOS; no libavahi headers needed.
#include <stdbool.h>

extern bool Static_avahi;

struct service_tab {
  void *buffer;
  char *line_type;
  char *interface;
  char *protocol;
  char *name;
  char *type;
  char *domain;
  char *dns_name;
  char *address;
  char *port;
  char *txt;
};

int avahi_browse(struct service_tab *table,int tabsize,char const *service_name);
void avahi_free_service_table(struct service_tab *table,int tabsize);

int avahi_start(char const *service_name,char const *service_type,int service_port,char const *dns_name,int base_address,char const *description);
#define AVAHI_H 1
#endif
