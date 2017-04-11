#include "../src/readBio.cpp"

int main(){
    readBio("/bioinfo/users/hcliment/easyGWASCore/data/testing/scones/skat/network.txt");

    string wd = "/bioinfo/users/hcliment/easyGWASCore/data/testing/scones/skat/";
    string pedBasename = wd + "genotype";
    string phenoFile = wd + "phenotype.txt";
    string netFile = wd + "network.txt";

    readBio(pedBasename, phenoFile, netFile, 0, 0.05);
}