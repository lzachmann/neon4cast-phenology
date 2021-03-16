#!/usr/bin/env bash

IMAGE=${1:-cpu}
TAG=${2:-latest}

GROUP=apis-staging
PROJECT=greenwave

if [ $IMAGE == "cpu" ]; then
		PORT=8888
		RUNTIME=""
		CRED=$PROJECT
		DATADRIVE="$(pwd)/data-raw"
elif [ $IMAGE == "gpu" ]; then
		PORT=8999
		RUNTIME="--runtime=nvidia"
		CRED=$PROJECT
		DATADRIVE="/datadrive"
else
    echo "Error: argument must be one of 'cpu' or 'gpu', you passed '$IMAGE'."
    exit 1
fi
ACCESS_POINT=http://localhost:$PORT/

IMAGE_NAME=registry.gitlab.com/$GROUP/$PROJECT/$IMAGE:$TAG
# IMAGE_NAME=rocker/ml:4.0.4

# PROFILE=~/profile/.rstudio
# mkdir -p $PROFILE
# cp -r $PROFILE "$PROFILE-$PROJECT"

if [[ "$(docker images -q $IMAGE_NAME 2> /dev/null)" == "" ]]; then
  echo "Image not found, building from recipe...."
  ./$(find . -path \*build.sh) $IMAGE $TAG
fi

echo "$ACCESS_POINT (with usr and pwd '$CRED')"

# this only works locally! on gpu it should be /.rstudio
docker run $RUNTIME -d -it --rm \
		--name rstudio-$PROJECT-$IMAGE-$TAG \
		-v "$(pwd)":/home/$CRED \
		-v "$(pwd)/cache":/tmp/cache \
		-v $DATADRIVE:/datadrive \
		-w /home/$CRED \
		-e USER=$CRED \
		-e PASSWORD=$CRED \
		-p $PORT:8787 \
		$IMAGE_NAME
# -v "$PROFILE-$PROJECT":"/home/$PROJECT/.rstudio" \

rm -rf kitematic/ #.rstudio/ rstudio/ "$PROFILE-$PROJECT"
