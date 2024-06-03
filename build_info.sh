if [ $1 = hash ]; then
    echo "EMME_COMMIT_HASH=\"$(git rev-parse HEAD)\""
fi
if [ $1 = time ]; then
    echo "EMME_BUILD_DATE=\"$(date -I'sec')\""
fi
