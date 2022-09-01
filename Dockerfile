FROM ubuntu:latest
MAINTAINER Dave Roe

# apt stuff
RUN apt-get update \
  && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip subversion gzip openjdk-17-jdk openjdk-17-doc groovy wget \
     zlib1g-dev gnuplot \
     bzip2 libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  && apt-get clean

# install stuff
ENV TMPDIR /tmp
RUN mkdir -p /opt/bin && cd /opt/bin \
  && wget -qO- http://get.nextflow.io | bash \
  && sed -i s/"curl -fsSL"/"curl -fsSLk"/ /opt/bin/nextflow \
  && chmod 755 /opt/bin/nextflow \
  && wget https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz \
  && gunzip KMC3.2.1.linux.tar.gz && tar -xvf KMC3.2.1.linux.tar && rm -f KMC3.2.1.linux.tar \
  && mv bin/* . && rmdir bin \
  && mkdir -p /opt/jars \
  && cd /opt/jars \
  && wget http://www.apache.org/dist/commons/math/binaries/commons-math3-3.6.1-bin.tar.gz \
  && tar -zxvf commons-math3-3.6.1-bin.tar.gz \
  && rm -f /opt/jars/commons-math3-3.6.1-bin.tar.gz
  
# google guava
ENV JAVA_HOME /usr/lib/jvm/java-17-openjdk-amd64
RUN cd /opt/jars \
  && wget https://repo1.maven.org/maven2/com/google/guava/guava/31.1-jre/guava-31.1-jre.jar \
  && chmod 755 guava-31.1-jre.jar

# env vars
ENV NXF_OPTS "-Xms4G -Xmx50G"
ENV JAVA_OPTS "-Xms4G -Xmx50G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH
ENV PATH /opt/bin:/opt/kpi:$PATH
ENV CLASSPATH /opt/jars/guava-31.1-jre.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:/usr/share/groovy/lib/commons-cli.jar:$CLASSPATH

# kpi files
RUN cd /opt && git clone https://github.com/droeatumn/kpi.git
CMD ["/opt/kpi/main.nf"]
