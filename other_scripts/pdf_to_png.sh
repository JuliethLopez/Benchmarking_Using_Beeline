for file in /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated/*.pdf; do
    filename=$(basename -- "$file")
    filename="${filename%.*}"
    echo "Converting $file to outputs_graphics/Curated/"$filename".png"
    pdftoppm -png -r 300 "$file" > "outputs_graphics/Curated/"$filename".png"
done