$0 ~ "^process" {
    gsub("^process | {", "", $0);
    proc_name=$0
}
$0 ~ "[ \t]+conda" {
    gsub("^[ ]+|conda |\"|\'", "", $0);
    print proc_name, $0
}