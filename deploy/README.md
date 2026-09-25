# Deploying `rms_gene_explorer` as a systemd service

Target: **Ubuntu 24.04, amd64**, Docker CE or `docker.io`, systemd 255.

## 1. Prerequisites on the server

```bash
docker --version                      # /usr/bin/docker must exist
systemctl is-active docker.service    # should print "active"
```

The repo must be checked out and the **store must already be built**. The store is not in
git and is not baked into the image — build it once on the server (or rsync it across):

```bash
cd ~/github_repos/2026-09-25_RNAseq_RMS_targets
Rscript bin/01_build_store.R
Rscript bin/02_repartition_store.R
Rscript bin/03_build_indexes.R
Rscript bin/04_validate_store.R      # must print "Store validation PASSED"
du -sh store                         # ~1.5 GB
```

Building the store needs the upstream `source_dir` from `config.yml` to be reachable. If
it is not, copy `store/` from a machine where it is — the app never reads `source_dir`.

## 2. Build the image

On amd64 no `--platform` flag is needed (that is only for Apple Silicon):

```bash
cd ~/github_repos/2026-09-25_RNAseq_RMS_targets
docker build -t rms-gene-explorer:4.5.3 .
```

The image is local-only and not in any registry, which is why the unit has **no
`docker pull`** step. After changing the Dockerfile, rebuild and
`sudo systemctl restart rms-gene-explorer`.

## 3. Install the unit

```bash
sudo install -m 0644 deploy/rms-gene-explorer.env.example /etc/rms-gene-explorer.env
sudo nano /etc/rms-gene-explorer.env          # set APP_DIR to the real checkout path
sudo install -m 0644 deploy/rms-gene-explorer.service /etc/systemd/system/
sudo systemctl daemon-reload
sudo systemctl enable --now rms-gene-explorer.service
```

`APP_DIR` must be a **literal absolute path**. systemd does not run a shell over the
environment file, so `$HOME` and `~` are not expanded and will silently break the mount.

## 4. Verify

```bash
systemctl status rms-gene-explorer.service
journalctl -u rms-gene-explorer.service -f          # container stdout, via journald
curl -so /dev/null -w '%{http_code}\n' http://127.0.0.1:3838/   # expect 200
docker exec rms-gene-explorer sh -c 'touch /app/store/x' # expect "Read-only file system"
```

## 5. Reaching the app

`BIND_ADDR` defaults to `127.0.0.1`, so the app is reachable only from the server itself.
Connect over an SSH tunnel from your laptop:

```bash
ssh -N -L 3838:127.0.0.1:3838 YOUR_USER@REMOTE_HOST
```

then open <http://127.0.0.1:3838>.

**Do not change `BIND_ADDR` to `0.0.0.0` casually.** PCGP is controlled-access under dbGaP
**phs000720** and the OPBG cohort is pre-publication; binding all interfaces publishes
sample-level expression to everything that can route to the host. If you need multi-user
access, put it behind institutional SSO or a reverse proxy with authentication, and check
it against the DUA first.

## 6. Operating it

| Task | Command |
|---|---|
| Restart after a code change | `sudo systemctl restart rms-gene-explorer` |
| Restart after a Dockerfile change | `docker build -t rms-gene-explorer:4.5.3 . && sudo systemctl restart rms-gene-explorer` |
| Stop and disable | `sudo systemctl disable --now rms-gene-explorer` |
| Logs since boot | `journalctl -u rms-gene-explorer -b` |
| Clear a crash-loop lockout | `sudo systemctl reset-failed rms-gene-explorer` |

## Design notes

**The container runs in the foreground.** No `-d`. With a detached container `ExecStart`
returns immediately, systemd concludes the service died, and `Restart=always` turns a
healthy app into a restart loop. Foreground plus `Type=simple` is what makes supervision
work.

**`ExecStartPre=-/usr/bin/docker rm -f %N`** clears a container orphaned by a hard reboot
or OOM kill. The leading `-` makes failure non-fatal, which matters on first start when
there is nothing to remove.

**Code directory writable, store read-only.** bslib compiles its Bootstrap theme into
`<appdir>/app_cache`; shiny fixes that path, so it cannot be redirected by an environment
variable, and a `--tmpfs` over it fails because the mountpoint cannot be created under a
read-only parent. Mounting the code writable and the store `:ro` protects the asset that
actually matters. Because the unit runs as root, `app_cache/` ends up root-owned inside the
checkout; it is gitignored and harmless, but `sudo chown -R $USER app_cache` if it is in
your way.

**`StartLimitBurst` / `StartLimitIntervalSec` are in `[Unit]`, not `[Service]`.** systemd
255 silently ignores them under `[Service]` — `systemd-analyze verify` reports
`Unknown key name`, and the crash-loop guard does nothing. Verified with
`systemctl show -p StartLimitBurst`, which reports `5` and `5min`.

**No `ProtectSystem=` / `PrivateTmp=` hardening.** This unit exists to talk to
`/run/docker.sock`. The usual sandboxing directives either block that outright or provide
no meaningful isolation given the daemon already runs as root. Harden the Docker daemon
instead.

**An alternative worth considering:** for a long-lived service, copying the app code
*into* the image and bind-mounting only `store:ro` would be stricter — immutable code,
one read-only mount, no `app_cache` write-back into the checkout, nothing writable in the
container. The cost is a rebuild on every code change. The current layout was kept because
it matches the mount location requested and lets you edit code and restart.
