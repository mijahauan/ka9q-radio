# Installing ka9q-radio on macOS

The primary platform for ka9q-radio is Debian Linux (see [INSTALL.md](INSTALL.md)), but the full package -- including the `radiod` server itself -- also builds and runs on macOS. Recent macOS SDKs provide everything the multicast core needs (`struct ip_mreqn`, `MCAST_JOIN_GROUP`), audio output uses PortAudio over CoreAudio, mDNS service advertising and discovery use the system's built-in mDNSResponder through the native `dns_sd` API (no Avahi needed), and the USB front ends are driven directly through libusb (no udev needed).

Differences from Linux to be aware of:

- No `SCHED_FIFO` real-time scheduling; `radiod` falls back to a raised nice priority. On Apple Silicon this has been fine in practice.
- Services are managed by launchd instead of systemd; use the provided `radiod-launchd` helper.
- The funcube dongle driver looks up its I/Q audio device by ALSA name and is not currently usable on macOS.
- SDRplay requires the proprietary SDRplay API; HydraSDR and Fobos require driver libraries not in Homebrew or MacPorts (build from source).

## Prerequisites

Install either MacPorts or Homebrew. If both are installed, the build prefers MacPorts; override with `make MACPREFIX=...`.

MacPorts:

```
sudo port install fftw-3-single iniparser libogg libopus libsamplerate libusb ncurses portaudio
sudo port install airspy airspyhf hackrf rtl-sdr    # optional, per your hardware
```

Homebrew:

```
brew bundle --file=src/Brewfile
```

Driver plugins for SDR front ends are enabled automatically at build time if the corresponding library headers are installed; `rx888.so`, `sig_gen.so`, and `funcube.so` need only the core libraries. Check the build output to see which drivers were built, or force one on/off with e.g. `make ENABLE_AIRSPY=1`.

## Compiling and Installing

Same as Linux:

```
cd ka9q-radio
make -j
sudo make install
```

On macOS this installs binaries to **/usr/local/sbin** and **/usr/local/bin**, driver plugins to **/usr/local/lib/ka9q-radio**, support files to **/usr/local/share/ka9q-radio**, and creates **/etc/radio**, **/var/lib/ka9q-radio**, **/etc/fftw**, and **/var/log/ka9q-radio**. The Linux-only pieces (systemd units, udev rules, sysctl.d, logrotate.d, cron.d, sysusers.d) are skipped; a launchd plist template and the `radiod-launchd` helper are installed instead. No `radio` user is created; run `radiod` as your own user (interactively) or as root (launchd).

## Configuring and Running radiod

Configuration is identical to Linux: create **/etc/radio/radiod@*instance*.conf** starting from the examples in the **config** subdirectory. See [ka9q-radio.md](ka9q-radio.md) and its companions.

Run it interactively first:

```
radiod -v /etc/radio/radiod@myinstance.conf
```

Then verify from another terminal that `control` discovers the instance (via mDNS/Bonjour) and that `monitor` plays audio. You can also check the advertisement directly with the system tool:

```
dns-sd -B _ka9q-ctl._udp
```

### Running under launchd

The equivalent of `sudo systemctl enable --now radiod@myinstance` is:

```
sudo radiod-launchd myinstance
```

This instantiates the plist template as **/Library/LaunchDaemons/net.ka9q.radiod.myinstance.plist**, loads it, and starts it. The service starts at boot, restarts if it dies, and logs to **/var/log/ka9q-radio/radiod@myinstance.log**. To stop and remove it:

```
sudo radiod-launchd myinstance remove
```

### RX888 firmware loading

There is no udev on macOS to trigger `rx888_boot`, and it isn't needed: set the `firmware=` option in the `[rx888]` hardware section of your config and `radiod` will upload the Cypress FX3 firmware itself when it finds an unloaded device:

```
[rx888]
device = rx888
firmware = SDDC_FX3.img
```

The image is installed in **/usr/local/share/ka9q-radio**. After a cold plug the first `radiod` start loads the firmware and the device re-enumerates; this takes a few extra seconds.

## Optimizations

### FFTW wisdom

As on Linux, `radiod` performance depends heavily on precomputed FFTW wisdom, especially at RX888 sample rates. See [FFTW3.md](FFTW3.md); `fft-gen` works the same on macOS and writes to **/var/lib/ka9q-radio**.

### Socket buffers

Linux installs get **/etc/sysctl.d/98-sockbuf.conf** to enlarge kernel socket buffers. There is no sysctl.d on macOS; the equivalent knob is `kern.ipc.maxsockbuf`. If `radiod` logs socket buffer warnings at high sample rates, raise it (the setting does not survive a reboot; add it to a launchd job if you need it permanently):

```
sudo sysctl -w kern.ipc.maxsockbuf=16777216
```

## Multicast notes

All ka9q-radio I/O is IP multicast, which works on macOS including over the loopback interface (`lo0` is multicast-capable out of the box, so the Linux `set_lo_multicast` service isn't needed for TTL=0 configurations). If you run `radiod` on a Mac with both Wi-Fi and Ethernet active, consider setting `iface=` in the `[global]` section to pin the output interface; Wi-Fi access points are generally poor at multicast.
