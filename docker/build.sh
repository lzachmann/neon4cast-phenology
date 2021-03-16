#!/usr/bin/env bash

IMAGE=${1:-cpu}
TAG=${2:-latest}

GROUP=apis-staging
PROJECT=$(basename $(pwd))

IMAGE_NAME=registry.gitlab.com/$GROUP/$PROJECT/$IMAGE:$TAG

FILE=$(find . -path \*$IMAGE/$TAG/Dockerfile)
docker build --no-cache --rm -t $IMAGE_NAME - < $FILE
