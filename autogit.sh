#! /bin/bash
if [ $# -ne 2 ];then echo -e "You need to provide comment & branch name !\n   autobash.sh \"MY COMMENT\" MYBRANCH\n" ; exit 0;fi

COMMENT="$1"
BRANCH=$2

git add . 
git commit -m"$COMMENT"
git push Nozzle $BRANCH
