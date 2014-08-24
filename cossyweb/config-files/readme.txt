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


############# Rserve configuration ####################
Rserv.conf: Rserve config file. location: /etc/Rserv.conf

## start RServe 
R CMD Rserve --RS-port 6312 --RS-workdir /etc/icossy/sessions/

## Stop RServe
Find the process id of Rserve: netstat -ap
Kill the process: sudo kill PID


############# reinstall icossy ########################
# detach icossy if it is being used
detach("package:icossy", unload=TRUE)

# uninstall
remove.packages("icossy")

# install
install.packages("PKG_FILE_PATH", repos = NULL, type = "source")
