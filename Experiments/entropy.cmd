openssl enc -aes-256-ctr -pass pass:"42" -nosalt </dev/zero 2>/dev/null | head -c 1024k > Experiments/entropy.out
