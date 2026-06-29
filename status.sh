#!/usr/bin/env bash
# Show status of the Coverage Checker web app.
set -uo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"
source "$DIR/server.conf"

PIDFILE="$DIR/server.pid"
LOG="$DIR/server.log"

echo "Coverage Checker — status"
echo "  Port:    $PORT"
echo "  Dir:     $DIR"

running=0
PID=""
if [[ -f "$PIDFILE" ]]; then
  PID="$(cat "$PIDFILE")"
  if kill -0 "$PID" 2>/dev/null; then
    running=1
  fi
fi

if [[ "$running" -eq 1 ]]; then
  echo "  State:   RUNNING (PID $PID)"
  # resource usage
  if command -v ps >/dev/null 2>&1; then
    read -r etime rss <<<"$(ps -o etime=,rss= -p "$PID" 2>/dev/null | awk '{print $1, $2}')"
    [[ -n "${etime:-}" ]] && echo "  Uptime:  $etime"
    [[ -n "${rss:-}" ]] && echo "  Memory:  $(( rss / 1024 )) MB"
  fi
else
  echo "  State:   STOPPED"
fi

# Health check (independent of PID file)
code=$(curl -s -o /dev/null -w "%{http_code}" "http://127.0.0.1:$PORT/" 2>/dev/null)
[[ -z "$code" ]] && code="000"
if [[ "$code" == "200" ]]; then
  echo "  Health:  OK (HTTP 200)"
  # how many panels are loaded
  n=$(curl -fsS "http://127.0.0.1:$PORT/api/panels" 2>/dev/null \
        | grep -o '"name"' | wc -l | tr -d ' ')
  [[ -n "$n" && "$n" != "0" ]] && echo "  Panels:  $n BED file(s) available"
  IP=$(hostname -I 2>/dev/null | awk '{print $1}')
  echo "  Local:   http://localhost:$PORT"
  [[ -n "${IP:-}" ]] && echo "  Network: http://$IP:$PORT"
else
  echo "  Health:  NOT RESPONDING (HTTP $code)"
fi

if [[ -f "$LOG" ]]; then
  echo "  Log:     $LOG"
  echo "  --- last 3 log lines ---"
  tail -n 3 "$LOG" | sed 's/^/  /'
fi

# Exit non-zero if not healthy (useful for monitoring)
[[ "$code" == "200" ]] && exit 0 || exit 1
