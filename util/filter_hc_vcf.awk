BEGIN {
    FS="\t"
    OFS="\t"
}

/^#/ {
    next
}

(($4 == "G" || $4 == "C" || $4 == "A" || $4 == "T") && 
 ($5 == "G" || $5 == "C" || $5 == "A" || $5 == "T")) {
    split($9, fields, ":")
    split($10, values, ":")
    for (i in fields)
        data[fields[i]] = values[i]

    if (data["GT"] != "0/1" && data["GT"] != "1/0")
        next

    split(data["AD"], ads, ",")
    print $1, $2, ads[1] + ads[2], ads[2]
}
