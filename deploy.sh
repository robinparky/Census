#alias lcbu='ant -f $HOME/pms/Census/build_simple.xml census'


ant -f  build_simple.xml  census
h=$(hostname)
if [[ $h == "yateslab100792.scripps.edu" ]]; then
	cp /home/yateslab/IdeaProjects/Census2/deploy/census2.jar /home/rpark/ip2_tomcat/webapps/ip2/softwares/census
	
fi
cp /home/rpark/git/Census/deploy/census2.jar /home/rpark/ip2_tomcat/webapps/ip2/softwares/census
#cp /home/rpark/git/Census/deploy/census2.jar /home/rpark/ip2_tomcat/webapps/ip2/WEB-INF/lib/census.jar
#cp /home/rpark/pms/Census/deploy/census2.jar /home/rpark/ip2_tomcat/webapps/ip2/WEB-INF/lib/census.jar



#cp -rf /home/rpark/pms/Census/build/classes/* /data/1/root/java/Census/.
#scp /home/rpark/pms/Census/deploy/census.jar root@shamu:/data/1/root/java/lib/census2_1.jar
#cp /home/rpark/pms/Census/deploy/census.jar /data/1/rpark/deploy/census.jar
#cp /home/rpark/pms/Census/deploy/census.jar /home/rpark/ip2_tomcat/webapps/ip2/softwares/census/census.jar
#scp /home/rpark/pms/Census/deploy/census.jar ip2@goldfish:/home/ip2/ip2_tomcat/webapps/ip2/softwares/census/.
