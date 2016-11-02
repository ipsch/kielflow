#!/bin/bash

#echo "OPTIND starts at $OPTIND"

# Global vars
filename=""
filename_out=""
target_dir="./"
script=""
output_local=false

function find_string () # with two different return values
# $1 - search string
# $2 - output (optional);
#      if not defined output to standard output
# Example for usage:
#   1) find_string '_data_splot_my' default_data
#   2) default_data=`find_string "no_string_correct"`
# IMPORTANT: never add any output / echo line to this
#            function if using option 2)
# Author: I.P.Schnell @ CAU Kiel (c) 2016
{
  local return_if="./data/splot_data.dat"
  local return_else="./data/plot_data.dat"
  local return_val=""
  
  flag=`echo $1|awk '{print match($1,"splot")}'`;
  #flag contains position of first appeareance of $1
  if [ "$flag" -gt 0 ];then
    return_val=$return_if
  else
    return_val=$return_else
  fi
  
  local __return__=$2 
  if [[ "$__return__" ]]; then # return_val to $2 if defined
    eval $__return__="'$return_val'"
  else # return_val to standard output if $2 not defined
    echo "$return_val"
  fi 
}



function plot_directive ()
# $1 - script
# $2 - data
# $3 - output

{
  local default_script="./plots/default_"$1".pl4"
  local custom_script="./plots/"$1".pl4"
  local custom_data=$2
  local default_data=""
  local plot_output=$3
  
  find_string $1 default_data
  
  rm $custom_script
  sed -e 's:'"$default_data"':'"$custom_data"':' $default_script >> $custom_script
  

  echo "my_plots: invoking plot4:"
  echo "  plot4 $custom_script -o $plot_output"
  plot4 "$custom_script" -o "$3"


}


######## MAIN ########

while getopts ":d:hp:" optname
  do
    case "$optname" in
    

    

 
      "d")
        target_dir=$OPTARG
        ;;
        

    
      "p")
        case "$OPTARG" in
          "plot_Ux") script="plot_Ux" ;;  
          "plot_ni") script="plot_ni" ;;
          "plot_Ph") script="plot_Ph" ;;
          "plot_force") script="plot_force" ;;  
          # 2d data
          "splot_ni") script="splot_ni" ;;
          "splot_Ph") script="splot_Ph" ;;
          "splot_FUx") script="splot_FUx" ;; 
          "splot_Fni") script="splot_Fni" ;;
          "splot_FPh") script="splot_FPh" ;;
          "splot_velocity") script="splot_velocity" ;;
          "splot_force") script="splot_force" ;;
          "?")
            echo "unshiny"
            ;;
        esac
        ;;
        
      "h")
        echo "-o | outpit file"
        echo "-p | plot Ux/Uy/Uz/ni/Ph"
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

shift $(($OPTIND - 1))

for file in "$@"; do  
  echo "my_plots: processing: $file"

  # get filename from input
  # dump directory
  t_file_=$(basename "$file")
  t_dir_=$(dirname "$file")
  
  # remove splot_data / plot_data substring
  # from filename if exists
  t_file_=${t_file_#splot_data *}
  t_file_=${t_file_#plot_data *}
  
  # add target directory (if specified)
  # append used scriptname to filename
  # change file extension from whatever extension to pdf
  target="$target_dir$script ${t_file_%.*}.pdf"
  
  echo "my_plots: target: $target"
  echo "my_plots: script: $script"

  plot_directive "$script" "$file" "$target"
  
done

  
  