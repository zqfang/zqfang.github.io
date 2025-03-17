# Cheatsheet for command line


usefull tools for linux command line

### Make terminal cool, install `OhMyZsh`

```shell
# install zsh
sudo apt-get install zsh # ubuntu
# change default shell to zsh
sudo chsh -s $(which zsh) $USER
## chsh -s /usr/bin/zsh 
# install ohmyzsh
sh -c "$(wget https://raw.github.com/robbyrussell/oh-my-zsh/master/tools/install.sh -O -)"

source ~/.zshrc
```

### ssh login withoutpassword

1. on your local machine, run
```bash
ssh-keygen
```
press Enter 3 time until finished. a `"~/.ssh/id_rsa.pub"` and `"~/.ssh.id_rsa"` will be generated.

2. copy the public key file to the remote machine
```bash
ssh-copy-id remote_username@remote_server_ip_address
```

3. then you can login without password
```bash
ssh remote_username@remote_server_ip_address
```


### Terminal keyboard shortcuts

- Jump to head: `Ctrl + a`  
- Jump to end: `Ctrl + e`  
- Delete strings ahead: `Ctrl + u`
- Delete strings follow: `Ctrl + k`

### Program keeps running in the background  

#### Run cmd using `nohup`

```shell
nohup command [options] &
```

#### Run cmd using `zellij`. 

`zellij` is drop-in replacement of tmux
```shell

# create a new session
zellij 
# specify a new session, named bio
zellij -s bio

# subcommds 
zellij ls # list all session
zellij a bio # attach session
zellij ka # kill all sessions
zellij k bio # kill a specfic session, bio
zellij da # delete all session
zellji d bio # detelete a specifc session, bio
```




#### Run cmd using `Tmux`

**Outside Tmux:**  
Typically these are run outside, but you can also run them inside an existing session  

a. Start New Session 

```shell
tmux new -s myname
```

b. Attach To Existing Session 

```shell
tmux attach -t myname          #by name 
tmux attach 4                  #by number (in this case 4)
```

c. List Sessions 
```shell
tmux ls
```

d. Kill Session
```shell 
tmux kill-session -t myname
```

**Inside Tmux Session:**  
Start each command with `CTRL + b`, release, then press one of the following:  

| Panes                      |                                                          |
| -------------------------  | -------------------------------------------------------- | 
| %                          | vertical split                                           |
| "                          | horizontal split                                         |
| d                          | detach from session (it keeps running in the background) |
| x                          | kill pane                                                |
| Up/Down/Left/Right         | move between panes                                       |
| PageUP/PageDown            | `CTRL+c` to exit the PageUp/Down mode                    |
| `Fn`+Up/Down               | PageUp/Down: Mac keyboard                                |
| : + resize-pane -D         | Resizes the current pane down                            |
| : + resize-pane -U         | Resizes the current pane upward                          |
| : + resize-pane -L         | Resizes the current pane left                            |
| : + resize-pane -R         | Resizes the current pane right                           |
| : + resize-pane -D 20      | Resizes the current pane down by 20 cells                |


### File compression and decompression

Decompression  

| File type         | Cmd                 | e.g.                              |
| ----------------- | ------------------- | --------------------------------- |
| *.tar             | tar -xvf            |                                   |
| *.tar.gz or *.tgz | tar -xvzf           |                                   |
| *bz2              | bzip2 -d or bunzip2 |                                   |
| *.tar.bz2         | tar -xjf            |                                   |
| *.Z               | uncompress          |                                   |
| *.tar.Z           | tar -xZf            |                                   |
| *.rar             | unrar e or rar x    | unrar e file.rar                  |
| *.zip             | unzip               |                                   |
| *.gz              | gunzip              |                                   |


Compression  

| File type         | Cmd                 | e.g.                              |
| ----------------- | ------------------- | --------------------------------- |
| *.tar             | tar -cvf            |                                   |
| *.tar.gz or *.tgz | tar -cvzf           |                                   |
| *bz2              | bzip2 -z            |                                   |
| *.tar.bz2         | tar -cjf            |                                   |
| *.Z               | compress            |                                   |
| *.tar.Z           | tar -cZf            |                                   |
| *.rar             | rar a               | rar a -ep1 newname /home/user/cpp |
| *.zip             | zip                 |                                   |
| *.gz              | gzip                |                                   |

For rar installation

```shell
sudo apt-get install rar
```

### Handy tricks for handling filepath

very useful to strip file sufix, path et.al.

```shell
# e.g.
var=./home/fastq/filename_R1.fq.gz

# extract filename
${var#*/} # -> home/fastq/filename_R1.fq.gz
var1=${var##*/}  # -> filename_R1.fq.gz

# remove file suffix
${var1%.*} # -> filename_R1.fq
${var1%%.*} # -> filename_R1

# get basebame
var2=$(basename "${var}" .fq.gz) #-> filename_R1
```


###  `vim` shortcuts

In view mode, curso movement
```
^ (shift + 6): jump to the start of current line
$ (shift + 4): jump to the end of current line

gg: go to the first line of the document
G: go to the last line of the document

H: move to the top of screen
M: move to the middle of screen
L: move to the bottom of screen
```

### 
