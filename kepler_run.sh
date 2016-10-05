#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH --time=6-00:00:00
#SBATCH --mail-type=END,FAIL                          # notifications for job done & fail
#SBATCH --mail-user=schnell@theo-physik.uni-kiel.de   # send-to address

default_file="./data/fields.h5"
argline=""

while getopts ":T:M:Q:" optname
  do
    case "$optname" in
      "T")
        argline="$argline -T $OPTARG"
        ;;
      "M")
        argline="$argline -M $OPTARG"
        ;;
      "Q")
        argline="$argline -Q $OPTARG"
        ;;
      "?")
        echo "Unknown option $OPTARG"
        ;;
        
      ":")
        echo "No argument value for option scale$OPTARG"
        ;;
      *)
      
      # Should not occur
        echo "Unknown error while processing options"
        ;;
    esac
done    

echo "argline:"
echo "$argline"

make clean
# input data
if [ -f ./data/fields.h5 ] ;
then
  echo "data was found at: $default_file"
  echo "using this as input data"
else
  echo "data was no data at: $default_file"
  echo "creating input using ./bin/frontend"
  make bin/frontend
./bin/frontend $argline
fi

# program execution
make clean
make bin/kielflow
./bin/kielflow $argline