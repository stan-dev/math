#!/usr/bin/python

# Script for retrieving release notes from pull requests merged since the last release

import urllib.request, json

release_url = "https://api.github.com/repos/stan-dev/math/releases/latest"
releases = json.loads(urllib.request.urlopen(release_url).read().decode())
last_release_date = releases["published_at"]

prs_url = "https://api.github.com/repos/stan-dev/math/pulls?state=closed&sort=created&per_page=100&sort=asc&page="


num_of_prs = 0
prs_on_page = 1
current_page = 1
# cycle through the pages until you hit a page with no PRs
while prs_on_page > 0:
    tmp_url = prs_url + str(current_page)
    with urllib.request.urlopen(tmp_url) as url:
        prs_info = json.loads(url.read().decode())
        prs_on_page = len(prs_info)
        for pr in prs_info:
            # check if PR was merged
            if pr["merged_at"]:
                # if merged check date
                if pr["merged_at"] > last_release_date:
                    num_of_prs = num_of_prs + 1
                    body = pr["body"]
                    parsing_release_notes = False
                    for line in body.split("\r\n"):
                        if parsing_release_notes:
                            if line.find("##") >= 0:
                                parsing_release_notes = False
                        if parsing_release_notes:
                            line = line.strip()
                            if len(line) > 0:
                                print(" - " + line + "(#" + str(pr["number"]) + ")")
                        if line.find("Release notes") >= 0:
                            parsing_release_notes = True
    current_page = current_page + 1

print("\nNumber of merged PRs:" + str(num_of_prs))
