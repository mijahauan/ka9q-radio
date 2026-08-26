# -*- makefile -*-
# sub-make for ka9q-radio systemd .service files
# Copyright 2026, Phil Karn, KA9Q

systemdunitdir ?= /etc/systemd/system
bindir         ?= /usr/local/bin
sbindir        ?= /usr/local/sbin
recordingsdir  ?= /srv/ka9q-radio/recordings
pkgdatadir     ?= /usr/local/share/ka9q-radio

UNAME_S := $(shell uname -s)

SERVICEFILES=acars.service aprsfeed.service cwd.service ft4-decode.service \
	ft8-decode.service ft4-decode@.service ft8-decode@.service \
	ft4-record.service ft8-record.service hfdl.service horusdemod.service \
	ka9q-radio@.service packetd.service radiod@.service recordings@.service \
	rx888_boot.service set_lo_multicast.service wspr-decoded.service

.PHONY: clean all install uninstall purge

ifeq ($(UNAME_S),Darwin)
# macOS has launchd, not systemd. Install the plist template (with install
# paths filled in); the radiod-launchd script instantiates it per instance.
all:
	@true

clean:
	@true

purge:
	@test "$(CONFIRMED)" = "1" || { \
	    echo "Refusing to run purge directly. Use top-level make purge."; \
	    exit 1; \
	}
	rm -f $(DESTDIR)$(pkgdatadir)/net.ka9q.radiod.plist.in

install:
	install -d -m 0755 $(DESTDIR)$(pkgdatadir)
	sed -e 's|@sbindir@|$(sbindir)|g' net.ka9q.radiod.plist.in \
	    > $(DESTDIR)$(pkgdatadir)/net.ka9q.radiod.plist.in
	chmod 644 $(DESTDIR)$(pkgdatadir)/net.ka9q.radiod.plist.in

uninstall:
	rm -f $(DESTDIR)$(pkgdatadir)/net.ka9q.radiod.plist.in

else
all: $(SERVICEFILES)

clean:
	rm -f *.service

purge:
	@test "$(CONFIRMED)" = "1" || { \
	    echo "Refusing to run purge directly. Use top-level make purge."; \
	    exit 1; \
	}
	if [ -n "$(DESTDIR)$(systemunitdir)" ]; then \
		cd $(DESTDIR)$(systemunitdir) && rm -f $(SERVICEFILES); \
	fi

install: $(SERVICEFILES)
	install -d -m 0755 $(DESTDIR)$(systemdunitdir)
	install -m 0644 $(SERVICEFILES) $(DESTDIR)$(systemdunitdir)
	if [ -z "$(DESTDIR)" ] && command -v systemctl >/dev/null; then \
		systemctl daemon-reload; \
	fi

uninstall:
	rm -f $(addprefix $(DESTDIR)$(systemdunitdir)/,$(SERVICEFILES))
	if [ -z "$(DESTDIR)" ] && command -v systemctl >/dev/null; then \
		systemctl daemon-reload; \
	fi
endif

%.service: %.service.in
	sed -e 's|@bindir@|$(bindir)|g' \
	    -e 's|@sbindir@|$(sbindir)|g' \
	    -e 's|@recordings@|$(recordingsdir)|g' \
	    $< > $@

show-service-paths:
	@echo bindir=$(bindir)
	@echo sbindir=$(sbindir)
	@echo systemdunitdir=$(systemdunitdir)
