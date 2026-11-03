#!/usr/bin/env bash
#
# 13_publish_datasette.sh — serve a built HPDB SQLite file through Datasette.
#
# The default mode starts a local read-only web interface so the seven tables
# can be browsed and queried in a browser. Set PUBLISH_TARGET to deploy a public
# instance instead (cloudrun, vercel or fly); those targets need the matching
# CLI and account credentials documented in the Datasette publishing guide.
#
#   bash server/13_publish_datasette.sh path/to/hpdb.sqlite
#   PUBLISH_TARGET=cloudrun SERVICE=hpdb bash server/13_publish_datasette.sh hpdb.sqlite
#
set -euo pipefail

DB_FILE="${1:?usage: $0 <hpdb.sqlite>}"
[ -f "$DB_FILE" ] || { echo "database not found: $DB_FILE" >&2; exit 1; }

HOST="${DATASSETTE_HOST:-127.0.0.1}"
PORT="${DATASSETTE_PORT:-8001}"
SERVICE="${SERVICE:-hpdb}"
PUBLISH_TARGET="${PUBLISH_TARGET:-local}"

if ! command -v datasette >/dev/null 2>&1; then
    echo "datasette not found. Install it with: pip install 'pufprobe[db]'" >&2
    exit 1
fi

# A short metadata block gives the instance a title and hides write access.
META_DIR="$(mktemp -d)"
META_FILE="$META_DIR/metadata.json"
cat > "$META_FILE" <<JSON
{
  "title": "Hypothetical Protein Database (HPDB)",
  "settings": {"allow_sql": true, "default_page_size": 50}
}
JSON
trap 'rm -rf "$META_DIR"' EXIT

case "$PUBLISH_TARGET" in
    local)
        echo "Serving $DB_FILE at http://$HOST:$PORT (Ctrl-C to stop)."
        exec datasette serve \
            --host "$HOST" --port "$PORT" \
            --metadata "$META_FILE" --immutable "$DB_FILE"
        ;;
    cloudrun|vercel|fly)
        echo "Publishing $DB_FILE to $PUBLISH_TARGET as $SERVICE."
        exec datasette publish "$PUBLISH_TARGET" "$DB_FILE" \
            --service "$SERVICE" --metadata "$META_FILE"
        ;;
    *)
        echo "unknown PUBLISH_TARGET '$PUBLISH_TARGET'; use local, cloudrun, vercel or fly." >&2
        exit 1
        ;;
esac
