add_block () {
    local name=$1
    cat << EOF > blocks/${name}
<+DESCRIPTION+>
<+UNIQUE/NONUNIQUE+>
<+DEPENDENCIES+>
<+OPTION+>
<+REQUIRED/OPTIONAL+>
<+TYPE+>
<+DESCRIPTION+>
<+REMARK+>
<+REMARK+>

<+OPTION+>
<+REQUIRED/OPTIONAL+>
<+TYPE+>
<+DESCRIPTION+>
<+REMARK+>
<+REMARK+>
EOF
}

add_option () {
    cat << EOF
<+OPTION+>
<+REQUIRED/OPTIONAL+>
<+TYPE+>
<+DESCRIPTION+>
<+REMARK+>
<+REMARK+>
EOF
}

edit_block () {
    vim blocks/$1
}

_cedit_block () {
    COMPREPLY=()
    for x in blocks/$2*
    do
        if [ -e $x ]
        then
            COMPREPLY+=("$(basename $x)")
        fi
    done
}

remove_block () {
    rm -i blocks/$1
}

complete -F _cedit_block edit_block
complete -F _cedit_block remove_block


