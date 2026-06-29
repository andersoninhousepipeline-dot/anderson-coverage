#!/usr/bin/env bash
# Stop the Coverage Checker web app.
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
source "$DIR/server.conf"

PIDFILE="$DIR/server.pid"
stopped=0

# 1) Stop via PID file
if [[ -f "$PIDFILE" ]]; then
  PID="$(cat "$PIDFILE")"
  if kill -0 "$PID" 2>/dev/null; then
    echo "Stopping PID $PID ..."
    kill "$PID" 2>/dev/null || true
    for i in $(seq 1 10); do
      kill -0 "$PID" 2>/dev/null || break
      sleep 0.5
    done
    if kill -0 "$PID" 2>/dev/null; then
      echo "  Did not exit; sending SIGKILL."
      kill -9 "$PID" 2>/dev/null || true
    fi
    stopped=1
  fi
  rm -f "$PIDFILE"
fi

# 2) Fallback: catch any stray app.py processes on our port
strays=$(pgrep -f "$PYTHON .*$APP" 2>/dev/null || true)
if [[ -n "$strays" ]]; then
  echo "Stopping stray process(es): $strays"
  kill $strays 2>/dev/null || true
  sleep 1
  kill -9 $strays 2>/dev/null || true
  stopped=1
fi

if [[ "$stopped" -eq 1 ]]; then
  echo "Stopped."
else
  echo "Not running."
fi
