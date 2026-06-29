#!/usr/bin/env bash
# Restart the Coverage Checker web app.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"$DIR/stop.sh"
sleep 1
exec "$DIR/start.sh"
