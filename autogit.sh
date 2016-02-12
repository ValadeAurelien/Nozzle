#! /bin/bash
exec git add . && echo "coucou"
exec git commit -m $1
exec git push Nozzle $2
