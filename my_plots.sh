#!/bin/bash
echo "OPTIND starts at $OPTIND"



filename=""

while getopts ":f:hp:" optname
  do
    case "$optname" in
    
      "f")
        echo "OPTARG is $OPTARG"
        filename=$OPTARG
        ;;
    
      "p")
        case "$OPTARG" in
          "Ux")
            konsole -e gnuplot -e "filename='$filename'" ./plots/interactive_Ux.pl
            ;;
          "Uy")
            echo "plot for Ux not implemented yet"
            echo "file: $filename"
            ;;
          "DC")
            rm ./plots/DC.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_DC.pl4 >> ./plots/DC.pl4
            plot4 ./plots/DC.pl4
            ;;
          "NC")
            rm ./plots/NC.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_NC.pl4 >> ./plots/NC.pl4
            plot4 ./plots/NC.pl4
            ;;
          "splot_density")
            rm ./plots/splot_density.pl4
            sed -e 's:"./data/splot_data.dat":"'$filename'":' ./plots/default_splot_density.pl4 >> ./plots/splot_density.pl4
            plot4 ./plots/splot_density.pl4
            ;;  
          "plot_density")
            rm ./plots/plot_density.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_plot_density.pl4 >> ./plots/plot_density.pl4
            plot4 ./plots/plot_density.pl4
            ;;
          "plot_potential")
            rm ./plots/plot_potential.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_plot_potential.pl4 >> ./plots/plot_potential.pl4
            plot4 ./plots/plot_potential.pl4
            ;;
          "splot_potential")
            rm ./plots/splot_potential.pl4
            sed -e 's:"./data/splot_data.dat":"'$filename'":' ./plots/default_splot_potential.pl4 >> ./plots/splot_potential.pl4
            plot4 ./plots/splot_potential.pl4
            ;;
          "plot_Ux")
            rm ./plots/plot_Ux.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_plot_Ux.pl4 >> ./plots/plot_Ux.pl4
            plot4 ./plots/plot_Ux.pl4
            ;;  
          "splot_velocity")
            rm ./plots/splot_velocity.pl4
            sed -e 's:"./data/splot_data.dat":"'$filename'":' ./plots/default_splot_velocity.pl4 >> ./plots/splot_velocity.pl4
            plot4 ./plots/splot_velocity.pl4
            ;;
          "plot_force")
            rm ./plots/plot_force.pl4
            sed -e 's:"./data/plot_data.dat":"'$filename'":' ./plots/default_plot_force.pl4 >> ./plots/plot_force.pl4
            plot4 ./plots/plot_force.pl4
            ;;  
          "splot_force")
            rm ./plots/splot_force.pl4
            sed -e 's:"./data/splot_data.dat":"'$filename'":' ./plots/default_splot_force.pl4 >> ./plots/splot_force.pl4
            plot4 ./plots/splot_force.pl4
            ;;
        esac
        ;;
    
        
        
        
      "h")
        echo "-f | filename"
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
    echo "OPTIND is now $OPTIND"
  done
  
  
  