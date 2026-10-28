"""Find the LMFDB label of a field we have computed.

A source is a function (degree, disc_abs) -> [(label, coeffs), ...] listing the LMFDB
fields of that degree and absolute discriminant. There are two:

  * api_source: the public LMFDB API (https://www.lmfdb.org/api/), no setup needed;
  * db_source:  a direct database connection, through the lmfdb package
                (`from lmfdb import db`) or lmfdb-lite (`from lmf import db`).
"""

import json
import time
import urllib.parse
import urllib.request

API_URL = "https://www.lmfdb.org/api/nf_fields/"
# The default Python-urllib user agent can be refused by the server.
USER_AGENT = "number-field-names (https://github.com/LMFDB/lmfdb/issues/6991)"


def api_source(degree, disc_abs, pause=0.2):
    rows = []
    offset = 0
    while True:
        query = urllib.parse.urlencode({
            "degree": f"i{degree}", "disc_abs": f"i{disc_abs}",
            "_fields": "label,coeffs", "_format": "json", "_offset": offset})
        url = f"{API_URL}?{query}"
        req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT,
                                                   "Accept": "application/json"})
        with urllib.request.urlopen(req) as resp:
            body = resp.read().decode("utf-8", "replace")
            ctype = resp.headers.get("Content-Type", "")
        try:
            data = json.loads(body)["data"]
        except (ValueError, KeyError):
            raise SystemExit(f"The LMFDB API did not return JSON data for\n  {url}\n"
                             f"Content-Type: {ctype}\nStart of the response:\n{body[:500]}\n\n"
                             "Try --source db instead (needs lmfdb-lite, see scripts/README.md).")
        rows += data
        time.sleep(pause)  # be gentle with the LMFDB server
        if len(data) < 100:  # the API returns at most 100 rows per page
            break
        offset += len(data)
    return [(r["label"], [int(c) for c in r["coeffs"]]) for r in rows]


def db_source():
    try:
        from lmfdb import db
    except ImportError:
        from lmf import db

    def source(degree, disc_abs):
        return [(r["label"], [int(c) for c in r["coeffs"]])
                for r in db.nf_fields.search({"degree": degree, "disc_abs": disc_abs},
                                             ["label", "coeffs"])]
    return source


def find_label(coeffs, r1, candidates, is_isomorphic):
    """Return the label among candidates [(label, coeffs)] of the field defined by coeffs
    with r1 real places, or None if it is not there.

    An exact match on the polredabs polynomial is tried first. In case PARI's polredabs
    differs from the one stored in the LMFDB, we fall back to an isomorphism test against
    the candidates with the same signature.
    """
    same_sig = [(label, c) for label, c in candidates if label.split(".")[1] == str(r1)]
    for label, c in same_sig:
        if tuple(c) == tuple(coeffs):
            return label
    for label, c in same_sig:
        if is_isomorphic(c, coeffs):
            return label
    return None
