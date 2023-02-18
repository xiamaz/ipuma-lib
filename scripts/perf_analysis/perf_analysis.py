import re
import json

from pathlib import Path

from argparse import ArgumentParser


def parse_logline(line):
  if "{" in line and "}" in line:
    start = line.find("{")
    end = line.rfind("}")
    entry = line[start:end+1]
    if m := re.match("] ([\w ]+) {"):
      tagname = m.group(1)
    else:
      tagname = ""
    return {
      "type": tagname,
      "data": json.loads(entry)
    }
  return None


def load_logfile(logfile_path):
  log_entries = []
  if logfile_path.exists():
    with open(logfile_path) as file:
      for line in file:
        if e := parse_logline(line):
          log_entries.append(e)
  return log_entries


parser = ArgumentParser()
parser.add_argument("logfile_path", type=Path)

args = parser.parse_args()

logs = load_logfile(args.logfile_path)