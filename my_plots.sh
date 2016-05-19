#!/bin/bash

#echo "OPTIND starts at $OPTIND"

# Global vars
filename=""
filename_out=""
script_name=""
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
  
  if [ "$filename_out" != "" ]; then
    echo "my_plots: invoking plot4:"
    echo "  plot4 $custom_script -o $plot_output"
    plot4 "$custom_script" -o "$plot_output"
  else
    echo "my_plots: invoking plot4:"
    echo "  plot4 $custom_script"
    plot4 "$custom_script"
  fi

}


######## MAIN ########

while getopts ":f:hlo:p:" optname
  do
    case "$optname" in
    
      "f")
        filename=$OPTARG
        ;;
    
      l)
	output_local=true
	;;
      
      o)
        filename_out=$OPTARG
        ;;
    
      "p")
        case "$OPTARG" in
          "plot_Ux") script_name="plot_Ux" ;;  
          "plot_density") script_name="plot_density" ;;
          "plot_potential") script_name="plot_potential" ;;
          "plot_force") script_name="plot_force" ;;  
          # 2d data
          "splot_density") script_name="splot_density" ;;
          "splot_potential") script_name="splot_potential" ;;
          "splot_velocity") script_name="splot_velocity" ;;
          "splot_force") script_name="splot_force" ;;
          "?")
            echo "unshiny"
            ;;
        esac
        ;;
        
      "h")
        echo "-f | filename"
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

# set output .pdf file to be located in same
# directory as .dat input files
if [ $output_local == true ]; then
  t_file_=$(basename "$filename")
  t_dir_=$(dirname "$filename")
  t_file_="${t_file_/.dat/.pdf}"
  t_file_="$script_name $t_file_"
  filename_out="$t_dir_/$t_file_"
fi


  #else

  #fi

echo "my_plots input: $filename"
echo "my_plots output: $filename_out"
echo "my_plots script: $script_name"

plot_directive "$script_name" "$filename" "$filename_out"

  
  