import re
import json


def parse_logline(line):
    if "{" in line and "}" in line:
        start = line.find("{")
        end = line.rfind("}")
        entry = line[start : end + 1]
        if m := re.search(r"\] ([\w ]+):? \{", line):
            tagname = m.group(1)
        else:
            tagname = ""
        try:
            entry_data = json.loads(entry)
            return {"type": tagname, "data": entry_data}
        except:
            print("Invalid line", line)


def multicmpHistogram(meta, match, line):
    if "Bucket" in line and "entries" in line:
        m = re.search(r"Bucket\s* (\d+)\s*has\s*(\d+)", line)
        if m is not None:
            bucket_id = m.group(1)
            entries = m.group(2)
            meta["histogram"][bucket_id] = entries


def totalExecutionTime(meta, match, line):
    if not "Local" in line:
        time_ms = match.group(1)
        meta["total_execution_time_ms"] = int(time_ms)


MATCHERS = {
    r"runIPUAlign@401": multicmpHistogram,
    r"^TIMEms:+\s*(\d+)": totalExecutionTime,
}


def match_meta(metadata, line):
    for regex, fun in MATCHERS.items():
        if m := re.search(regex, line):
            fun(metadata, m, line)


def load_logfile(logfile_path):
    log_entries = []
    if logfile_path.exists():
        metadata = {
            "name": logfile_path.stem,
            "histogram": {}
        }

        with open(logfile_path) as file:
            for line in file:
                if e := parse_logline(line):
                    log_entries.append(e)
                match_meta(metadata, line)
        metadata["entries"] = log_entries
        return metadata
