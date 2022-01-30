import os
from collections import defaultdict
from time import sleep
import sys
from getpass import getpass

if __name__ == "__main__":

    if len(sys.argv) > 1:

        gmt_files = sys.argv[1]
        files = []
        with open(gmt_files) as f:

            for line in f:

                files.append(line[:-1])
        files = sorted(files)

    else:

        files = sorted(["./gmtFiles/" + f for f in os.listdir("gmtFiles")])

    prefixes = {f[:-6] for f in files}
    pairs = defaultdict(list)
    to_wait = 60*45

    user_key = getpass("Enter user key:")

    for f in files:

        pairs[f[:-6]].append(f)

    pairs = [p for p in pairs.values() if len(p) == 2] # in case one set is missing
    
    for i, (DN,UP) in enumerate(pairs,start = 1):

        print(DN, UP)
        print(":")
        name = DN.split("/")[-1].split("_")[0]
        print(name, UP, DN)
        cmd = f'curl -i -X POST -H "user_key:{user_key}"  -H "Content-Type: multipart/form-data" -F "tool_id=sig_gutc_tool" -F "uptag-cmapfile=@{UP}" -F "dntag-cmapfile=@{DN}" -F "name={name}" -F "data_type=L1000" -F "dataset=Touchstone" https://api.clue.io/api/jobs > CMAP_responses/{name}_CMAP_response'
        os.system(cmd)
        if i % 10 == 0:

            sleep(to_wait)

