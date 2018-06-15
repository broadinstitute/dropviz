#!/bin/sh

# fuse mount point must happen during docker run, not during build, because there's no device or privs then
# despite docs, then doesn't work because it requires root, but gcsfuse isn't accessible by non-root processes if mounted with root
#mount /var/lib/shiny-server/bookmarks/

# must run as shiny https://github.com/GoogleCloudPlatform/gcsfuse/issues/175
gcsfuse -o allow_other -file-mode=777 -dir-mode=777 dropviz-bookmarks /var/lib/shiny-server/bookmarks

exec shiny-server > /var/log/shiny-server.log 2>&1

