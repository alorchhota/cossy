############## R Requirement ###########################
1) Install R (>=2.15.3)
2) Install packages: limma, igraph, jsonlite, Rserve, cossy, icossy

############# Directory requirement ####################
1) create a directory: /etc/icossy/sessions
2) copy all files and folders from ./cossyweb/php to /var/www/icossy/

############# apache2 configuration ####################
icossy: site config file. location: /etc/apache2/sites-available/icossy

## Enable/disable icossy site
sudo a2ensite icossy
sudo a2dissite icossy

## start, restart, stop, and status of apache2
sudo /etc/init.d/apache2 start
sudo /etc/init.d/apache2 restart
sudo /etc/init.d/apache2 stop
sudo /etc/init.d/apache2 status

## remember to add port number in /etc/apache2/ports.conf
Listen 8100

## remember to have sufficient upload-size
## edit php.ini file (possibly in /etc/php5/apache2/)
upload_max_filesize = 30M

############# Rserve configuration ####################
Rserv.conf: Rserve config file. location: /etc/Rserv.conf

## start RServe 
R CMD Rserve --RS-port 6312 --RS-workdir /etc/icossy/sessions/

## Stop RServe
Find the process id of Rserve: netstat -ap
Kill the process: sudo kill PID


## To access files in apache from Rserve
## or to have write permission in the Rserve folder
## both Rserve and apche has to be run with same user and group.

## Apache Config: add the following lines in /etc/apache2/envvars
export APACHE_RUN_USER=ashis
export APACHE_RUN_GROUP=ashis

## you may have to change the ownership of thevar/locks/apache2 folder.
sudo chown -R ashis:ashis /var/locks/apache2/

## Rserve config: Find the uid and gid of the user (here ashis:ashis)
id ashis

## then, edit the file /etc/Rserve.conf
uid UID_OF_USER
gid GID_OF_USER

## You have to restart both apache and Rserve.



############# reinstall icossy ########################
# detach icossy if it is being used
detach("package:icossy", unload=TRUE)

# uninstall
remove.packages("icossy")

# install
install.packages("PKG_FILE_PATH", repos = NULL, type = "source")
