#!/bin/bash

# configure your environment variable file found at minotourapp/envs.sh
source envs.sh
NAME=$1
WORKER_COUNT=$MT_CELERY_WORKER_COUNT

# If name was not specified use $USER                                           
[ -z "$NAME" ] && NAME=$USER

# If tmux has a session with name, attach that and exit                         
tmux has-session -t $NAME 2> /dev/null                                          
if [ $? -eq 0 ]; then                                                            
    tmux attach -t $NAME                                                        
    exit 0;                                                                     
fi

tmux new-session -d -s $NAME \;
tmux new-window -a -n dshell \;
tmux select-window -t dshell \;
#tmux send-keys "cd .." C-m \;
# If the python environment doesn't exist
tmux send-keys "source extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "python manage.py shell_plus" C-m \; 
tmux new-window -a -n minotour \; 
tmux select-window -t minotour \; 
#tmux send-keys "cd .." C-m \;
tmux send-keys "source extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "python manage.py runserver 8100" C-m \;
tmux split-window -v -p 70 \; 
#tmux send-keys "cd .." C-m \;
tmux send-keys "source extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "celery -A minotourapp worker -l info -f ~/data/logs/celery.log --concurrency $WORKER_COUNT -Ofair" C-m \;
tmux split-window -h \; 
tmux send-keys "mtlog" C-m \;
tmux new-window -a -n zshelly \;
tmux select-window -t zshelly\;
#tmux send-keys "cd .." C-m \;
tmux send-keys "source ../extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux new-window -a -n minimap2Celery \;
tmux select-window -t minimap2Celery \;
#tmux send-keys "cd .." C-m \;
tmux send-keys "source ../extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "celery -A minotourapp worker -l info -f $MT_LOG_FOLDER/celery.log -n minimap2 --concurrency $MT_CELERY_MINIMAP2_WORKER_COUNT -Ofair -Q minimap --pool threads" C-m \;
tmux split-window -h \;
#tmux send-keys "cd .." C-m \;
tmux send-keys "source ../extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "celery beat -A minotourapp -l info" C-m \;
tmux split-window -v \;
tmux send-keys "mtlog | grep align" C-m \;
tmux new-window -a -n mysql \;
tmux select-window -t mysql \;
tmux send-keys "mysql -u root -p issue_192" C-m \;
tmux new-window -a -n flower \;
tmux select-window -t flower \;
#tmux send-keys "cd .." C-m \;
tmux send-keys "source ../extra/minotour/bin/activate" C-m \;
tmux send-keys "./envs.sh" C-m \;
tmux send-keys "flower -A minotourapp --port=5556" C-m \;

tmux attach -t $NAME
