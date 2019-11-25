ant yguard 
cd deploy/zip
cp ../census_obf.jar census.jar 
java -jar census.jar -v
zip census_new *
#scp census_new.zip tomcat@bart.scripps.edu:/home/tomcat/tomcat/webapps/pmsws/project/census_data/data/deploy/.
scp census_new.zip tomcat@bart.scripps.edu:.
ssh tomcat@bart
