
for VARIABLE in "C" "cc" "hh" "h" "cpp"
do
    find . -name "*.${VARIABLE}" -print0 | while read -d $'\0' file
    do
            echo "${file}"
            clang-format -style="{BasedOnStyle: Google,SortIncludes: false}" -i ${file}
    done
done
