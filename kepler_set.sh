# upload: local -> remote



user='schnell'
server='Kepler'  # Kepler : 134.245.67.20

read -p "enter password @ "$server": " password

dir_base=$(date +"%y-%m-%d_%H%M%S")
dir_sub=$dir_base/{config,data}                       # subdirectory structure (example: $dir_base/{src,doc/{api,system},tools,db})



copy_source='*.cpp *.hpp makefile script_launch@kepler.sh'
copy_dir='config'




# create file structure
echo 'create file structure'

sshpass -p $password ssh "$server" "mkdir -p $dir_base" 
sshpass -p $password ssh "$server" "mkdir -p $dir_sub" 



# copy files and subdirectory content
echo 'copy files and subdirectory content'


sshpass -p $password scp -r $copy_dir    $user@$server:$dir_base
sshpass -p $password scp -r $copy_source $user@$server:$dir_base
sshpass -p $password scp ./data/fields.dat schnell@kepler:./$dir_base/data/fields.dat
sshpass -p $password ssh "$user@$server" "sbatch -D ./$dir_base ./$dir_base/script_launch@kepler.sh "


echo 'done'