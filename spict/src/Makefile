all: spict.so

spict.so: spict.cpp
	Rscript --vanilla -e 'TMB::compile("spict.cpp")'

clean:
	rm -rf *o
