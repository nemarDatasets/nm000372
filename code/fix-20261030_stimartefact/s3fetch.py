"""Fetch annexed files of a published NEMAR dataset by git-annex key from the public S3 URL.

materialize(repo, relpaths, dest) -> dict relpath -> local path. Each annexed file is downloaded from
https://nemar.s3.us-east-2.amazonaws.com/<id>/objects/<KEY> and its size + sha256 are checked against the key.
Plain git files are copied. Nothing is written to the repo.
"""
import hashlib, os, re, shutil, urllib.request

BASE = "https://nemar.s3.us-east-2.amazonaws.com"
KEYRE = re.compile(rb"SHA256E-s(\d+)--([0-9a-f]{64})[^\s/]*")


def annex_key(path):
    with open(path, "rb") as f:
        head = f.read(300)
    if head.startswith(b"/annex/objects/") or b"/annex/objects/" in head[:200]:
        m = KEYRE.search(head)
        return m.group(0).decode(), int(m.group(1)), m.group(2).decode()
    if os.path.islink(path):
        t = os.readlink(path).encode()
        m = KEYRE.search(t)
        if m:
            return m.group(0).decode(), int(m.group(1)), m.group(2).decode()
    return None


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for c in iter(lambda: f.read(1 << 24), b""):
            h.update(c)
    return h.hexdigest()


def materialize(repo, relpaths, dest, ds=None):
    ds = ds or os.path.basename(os.path.normpath(repo))
    out = {}
    for rel in relpaths:
        src = os.path.join(repo, rel)
        dst = os.path.join(dest, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        k = annex_key(src)
        if k is None:
            shutil.copyfile(src, dst)
        else:
            key, size, digest = k
            if not (os.path.exists(dst) and os.path.getsize(dst) == size and sha256(dst) == digest):
                url = f"{BASE}/{ds}/objects/{key}"
                tmp = dst + ".part"
                for attempt in range(4):
                    try:
                        with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
                            shutil.copyfileobj(r, f, 1 << 22)
                        break
                    except Exception as e:  # retry transient errors
                        if attempt == 3:
                            raise
                got = sha256(tmp)
                if os.path.getsize(tmp) != size or got != digest:
                    raise RuntimeError(f"{rel}: size/sha mismatch {os.path.getsize(tmp)} {got} vs {size} {digest}")
                os.replace(tmp, dst)
        out[rel] = dst
    return out
