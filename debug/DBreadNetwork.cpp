#include "../src/readGWAS.cpp"

int main(){

    std::string wd = "src/scones2/data/testing/scones/skat/";
    std::string pedBasename = wd + "genotype";
    std::string phenoFile = wd + "phenotype.txt";
    std::string netFile = wd + "network.txt";

    readGWAS(pedBasename, phenoFile, netFile, 0, 0.05);
}