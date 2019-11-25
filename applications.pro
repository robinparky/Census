#
# This ProGuard configuration file illustrates how to process applications.
# Usage:
#     java -jar proguard.jar @applications.pro
#

# Specify the input jars, output jars, and library jars.

-injars  deploy/censuscore.jar
-outjars deploy/censuscore_obf.jar

-libraryjars /usr/java/jdk/lib/dt.jar
-libraryjars lib/AbsoluteLayout.jar  
-libraryjars lib/commons-cli-1.0.jar
-libraryjars lib/jdom.jar
-libraryjars lib/junit-3.8.1.jar
-libraryjars lib/mzxmlviewer.jar
-libraryjars lib/swing-layout-0.7.jar
-libraryjars lib/TableLayout.jar
-libraryjars lib/trove.jar

#-libraryjars <java.home>/lib/rt.jar
#-libraryjars junit.jar
#-libraryjars servlet.jar
#-libraryjars jai_core.jar
#...

# Preserve all public applications.

-keepclasseswithmembers public class * {
    public static void main(java.lang.String[]);
}

# Print out a list of what we're preserving.

-printseeds

# Preserve all annotations.

-keep class * {
    public ** **;
}
