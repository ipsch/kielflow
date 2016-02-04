# download: remote -> local

user='schnell'
server='kepler'

read -p "enter password @ "$server": " password

sshpass -p $password scp -r $user@$server:./ ./data_kepler/
# sshpass -p $password ssh  "$server" "rm -r *"


