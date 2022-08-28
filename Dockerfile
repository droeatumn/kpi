FROM ubuntu:latest
MAINTAINER Dave Roe

# apt stuff
RUN apt-get update \
  && apt-get install -qyy curl git make vim cmake \
     gcc g++ unzip subversion gzip openjdk-8-jdk openjdk-8-doc groovy wget \
     zlib1g-dev gnuplot lynx maven \
     bzip2 libbz2-dev liblzma-dev libncurses5-dev libncursesw5-dev \
  && apt-get clean

# install stuff
ENV TMPDIR /tmp
RUN mkdir -p /opt/bin && cd /opt/bin \
  && wget -qO- http://get.nextflow.io | bash \
  && sed -i s/"curl -fsSL"/"curl -fsSLk"/ /opt/bin/nextflow \
  && chmod 755 /opt/bin/nextflow \
  && /opt/bin/nextflow \
  && wget https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz \
  && gunzip KMC3.2.1.linux.tar.gz && tar -xvf KMC3.2.1.linux.tar && rm -f KMC3.2.1.linux.tar \
  && mv bin/* . && rmdir bin \
  && mkdir -p /opt/jars \
  && cd /opt/jars \
  && wget http://www.apache.org/dist/commons/math/binaries/commons-math3-3.6.1-bin.tar.gz \
  && tar -zxvf commons-math3-3.6.1-bin.tar.gz \
  && rm -f /opt/jars/commons-math3-3.6.1-bin.tar.gz
  
# google guava
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64/
RUN cd /opt \
  && git clone https://github.com/google/guava.git \
  && cd guava/guava \
  && mvn install

# env vars
ENV NXF_OPTS "-Xms4G -Xmx50G"
ENV JAVA_OPTS "-Xms4G -Xmx50G"
ENV LD_LIBRARY_PATH /opt/lib:$LD_LIBRARY_PATH
ENV PATH /opt/bin:$PATH
ENV CLASSPATH /opt/guava/guava/target/guava-HEAD-jre-SNAPSHOT.jar:/opt/jars/commons-math3-3.6.1/commons-math3-3.6.1.jar:$CLASSPATH

# kpi files
RUN cd /opt && git clone https://github.com/droeatumn/kpi.git
CMD ["/opt/kpi/main.nf"]
