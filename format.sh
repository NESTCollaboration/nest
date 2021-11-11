FORMAT='clang-format -style="{BasedOnStyle: Google,SortIncludes: false}" -i'

for VARIABLE in "C" "cc" "hh" "h" "cpp"
do
    find . -type f -name "*.$VARIABLE" -exec echo ${FORMAT} \{\} \;
done

