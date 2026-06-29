#!/usr/bin/env bash
# Start the Coverage Checker web app (idempotent; writes a PID file).
set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
source "$DIR/server.conf"

PIDFILE="$DIR/server.pid"
LOG="$DIR/server.log"

# Already running?
if [[ -f "$PIDFILE" ]] && kill -0 "$(cat "$PIDFILE")" 2>/dev/null; then
  echo "Already running (PID $(cat "$PIDFILE")) on port $PORT."
  echo "  URL: http://localhost:$PORT"
  exit 0
fi
rm -f "$PIDFILE"

# Port already in use by something else?
if command -v ss >/dev/null 2>&1 && ss -ltn 2>/dev/null | grep -q ":$PORT[[:space:]]"; then
  echo "ERROR: port $PORT is already in use by another process." >&2
  echo "       Set a different port:  PORT=8101 ./start.sh" >&2
  exit 1
fi

echo "Starting Coverage Checker on $HOST:$PORT ..."
PORT="$PORT" nohup "$PYTHON" "$APP" >> "$LOG" 2>&1 &
PID=$!
echo "$PID" > "$PIDFILE"

# Wait for it to come up (panel index loads in ~3s)
for i in $(seq 1 40); do
  if ! kill -0 "$PID" 2>/dev/null; then
    echo "ERROR: process exited during startup. Last log lines:" >&2
    tail -n 8 "$LOG" >&2
    rm -f "$PIDFILE"
    exit 1
  fi
  if curl -fsS -o /dev/null "http://127.0.0.1:$PORT/" 2>/dev/null; then
    echo "Started (PID $PID)."
    IP=$(hostname -I 2>/dev/null | awk '{print $1}')
    echo "  Local:   http://localhost:$PORT"
    [[ -n "${IP:-}" ]] && echo "  Network: http://$IP:$PORT"
    echo "  Log:     $LOG"
    exit 0
  fi
  sleep 1
done

echo "WARNING: started (PID $PID) but did not pass health check in time." >&2
echo "Check the log: $LOG" >&2
exit 1
