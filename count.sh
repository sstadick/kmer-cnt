#!/usr/bin/sh

echo "Timing moded2 version"
/usr/bin/time bash -c "zcat M_abscessus_HiSeq_10M.fa.gz | head -n 1000000 | ./kc-py-mod2/target/release/kc-py-mod2 >/dev/null"

echo "Timing moded version"
/usr/bin/time bash -c "zcat M_abscessus_HiSeq_10M.fa.gz | head -n 1000000 | ./kc-py-mod/target/release/kc-py-mod >/dev/null"

echo "Timing single threaded"
/usr/bin/time bash -c "zcat M_abscessus_HiSeq_10M.fa.gz | head -n 1000000 | ./kc-py/target/release/kc-py >/dev/null"

