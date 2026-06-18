"""Read and write workspace assignments in Hermes Desktop LevelDB."""
import plyvel, json, sys

db_path = f"C:/Users/朱浩钒/AppData/Roaming/hermes/Local Storage/leveldb"
current_session = "20260614_141359_c970ae47"  # This session

try:
    db = plyvel.DB(db_path, create_if_missing=False)
except Exception as e:
    print(f"Error opening DB: {e}")
    sys.exit(1)

# Search for workspace keys
for key_bytes, val_bytes in db:
    key = key_bytes.decode('utf-8', errors='replace')
    val = val_bytes.decode('utf-8', errors='replace')
    if 'workspace' in key.lower() or 'mnt' in key.lower():
        print(f"Key: {key}")
        print(f"Val: {val[:200]}")
        print("---")

db.close()
