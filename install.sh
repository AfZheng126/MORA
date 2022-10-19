mkdir binaries

# Installing Pufferfish

echo ""
echo "Building Pufferfish"
echo ""

cd binaries
git clone https://github.com/COMBINE-lab/pufferfish.git
cd pufferfish
mkdir build
cd build
cmake ../
make
cd src
cp pufferfish ../../
cd ../../../

# Installing Bowtie2
echo ""
echo Building Bowtie2
echo ""

git clone https://github.com/BenLangmead/bowtie2.git
cd bowtie2
make
cd ../

# Installing Minimap2
echo ""
echo Building Minimap2
echo ""

git clone https://github.com/lh3/minimap2
cd minimap2
make
cd ../../


#Give executable permission to scripts
cd binaries
chmod +x pufferfish/pufferfish
chmod +x bowtie2/bowtie2
chmod +x bowtie2/bowtie2-build
chmod +x minimap2/minimap2
