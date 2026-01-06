use crate::cli::ClusterArgs as Cli;
use fastbloom::BloomFilter;
use std::sync::Arc;
use std::sync::Mutex;
use std::thread;
use fxhash::FxHashMap;
use crate::seeding;
use crate::utils::*;
use std::io::BufReader;
use crossbeam_channel::bounded;

pub fn read_to_split_kmers(
    k: usize,
    l: usize,
    threads: usize,
    args: &Cli,
) -> (Vec<(u64,[u32;2])>, Vec<(u64,[u32;2])>){


    let (bf_snpmer_maps, bf_blockmer_maps) = first_iteration(k, l, threads, args);
    //log::info!("Finished with bloom filter processing in {:?}. Round 2 - Start k-mer counting...", start.elapsed());
    //log_memory_usage(true, "Memory usage after bloom filter processing");

    let start = std::time::Instant::now();
    let (snpmer_maps, blockmer_maps) = second_iteration(k, l, threads, args, bf_snpmer_maps, bf_blockmer_maps);
    let snpmer_size_raw = snpmer_maps.iter().map(|x| x.len()).sum::<usize>();
    log::info!("Finished with k-mer counting in {:?}. Total SNPmers after bloom filter: {}",
        start.elapsed(), snpmer_size_raw);
    log_memory_usage(true, "Memory usage after second round of k-mer counting processing");

    // Vectorize and filter SNPmers
    let mut vectorized_snpmers = vec![];
    for map in snpmer_maps.into_iter(){
        for (kmer, counts) in map.into_iter(){
            if args.single_strand {
                if counts[0] > 2 {
                    vectorized_snpmers.push((kmer, counts));
                }
            }
            else{
                if counts[0] > 0 && counts[1] > 0 && counts[0] + counts[1] > 2{
                    vectorized_snpmers.push((kmer, counts));
                }
            }
        }
    }

    // Vectorize and filter blockmers
    let mut vectorized_blockmers = vec![];
    for map in blockmer_maps.into_iter(){
        for (kmer, counts) in map.into_iter(){
            if args.single_strand {
                if counts[0] > 2 {
                    vectorized_blockmers.push((kmer, counts));
                }
            }
            else{
                if counts[0] > 0 && counts[1] > 0 && counts[0] + counts[1] > 2{
                    vectorized_blockmers.push((kmer, counts));
                }
            }
        }
    }

    let snpmer_size_retain = vectorized_snpmers.len();
    let blockmer_size_retain = vectorized_blockmers.len();
    log::info!("Removed {} SNPmers with counts < 1 in both strands and <= 3 multiplicity.", snpmer_size_raw - snpmer_size_retain);
    //log::info!("Removed {} blockmers with counts < 1 in both strands and <= 3 multiplicity.", blockmer_size_raw - blockmer_size_retain);
    if snpmer_size_retain < snpmer_size_raw / 1000 {
        log::error!("Less than 0.1% of SNPmers have counts > 1 in both strands and > 2 multiplicity. This may indicate a problem with the input data or very low coverage. Consider using --single-strand");
        std::process::exit(1);
    }
    log::debug!("Final SNPmer hashmap len after vectorization: {}", snpmer_size_retain);
    log::debug!("Final blockmer hashmap len after vectorization: {}", blockmer_size_retain);
    log_memory_usage(false, "Memory usage after second round of k-mer counting processing");

    return (vectorized_snpmers, vectorized_blockmers);
}

fn first_iteration(
    k: usize,
    l: usize,
    threads: usize,
    args: &Cli,
) -> (Vec<FxHashMap<u64,[u32;2]>>, Vec<FxHashMap<u64,[u32;2]>>)
{
    //Topology is
    //      A-SEND: tx_head, , B-REC: rx_head1, rx_head2...
    // |   |  ...
    // B   B  ...  B-SEND: txs[0...], txs2[0...],... C-REC: rxs
    // | x | x | ...
    // C   C  ...
    let hm_size = threads;
    let mask = !(1 << 63);
    let bf_size = args.bloom_filter_size;
    let aggressive_bloom = args.aggressive_bloom;
    let mut bf_snpmer_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); hm_size];
    let mut bf_blockmer_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); hm_size];
    if bf_size > 0.{
        let num_b = threads/10 + 1;
        let counter = Arc::new(Mutex::new(0));
        let mut rxs = vec![];
        let mut txs_vecs = vec![vec![]; num_b];
        for _ in 0..threads {
            //let (tx, rx) = unbounded();
            let (tx, rx) = bounded(500);
            for i in 1..num_b{
                txs_vecs[i].push(tx.clone());
            }
            txs_vecs[0].push(tx);
            rxs.push(rx);
        }

        //let (tx_head, rx_head1) = unbounded();
        let (tx_head, rx_head1) = bounded(500);
        let mut rx_heads = vec![];
        for _ in 1..num_b{
            let rx_head2 = rx_head1.clone();
            rx_heads.push(rx_head2);
        }
        rx_heads.push(rx_head1);

        assert!(txs_vecs.len() == rx_heads.len());

        let fq_files = args.input_files.clone();
        let minimum_base_quality = args.minimum_base_quality;
        let use_blockmers = args.use_blockmers;
        //A: Get k-mers
        thread::spawn(move || {
            for fq_file in fq_files{
                let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
                let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
                while let Some(record) = reader.next() {
                    let rec = record.expect("Error reading record");
                    let rec_id_string = std::str::from_utf8(rec.id()).unwrap();
                    let last_field = rec_id_string.split_whitespace().last().unwrap_or(""); 
                    let seq;
                    let qualities; 
                    if last_field == "rc"{
                        seq = reverse_complement(&rec.seq().to_vec());
                        qualities = rec.qual().map(|q| q.iter().rev().cloned().collect());
                    }
                    else{
                        seq = rec.seq().to_vec();
                        qualities = rec.qual().map(Vec::from);
                    }
                    tx_head.send((seq,qualities)).unwrap();
                }
            }
            drop(tx_head);
            log::debug!("Finished reading all reads.");
        });
        //B: Process kmers and send to hash maps
        for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()){
            let clone_counter = Arc::clone(&counter);
            thread::spawn(move || {
                loop{
                    match rx_head.recv() {
                        Ok((seq, qualities)) => {
                            // Extract both SNPmers and blockmers
                            let split_kmer_info = seeding::split_kmer_mid(seq.clone(), qualities.clone(), k, minimum_base_quality);
                            let blockmer_info = seeding::blockmer_kmers(seq, qualities, k, l, minimum_base_quality);

                            // Hash and distribute SNPmers
                            let mut snpmer_vec_and_canon = vec![vec![]; hm_size];
                            for kmer_i_and_canon in split_kmer_info.into_iter() {
                                let kmer = kmer_i_and_canon & mask;
                                let hash = kmer % threads as u64;
                                snpmer_vec_and_canon[hash as usize].push(kmer_i_and_canon);
                            }

                            // Hash and distribute blockmers
                            let mut blockmer_vec_and_canon = vec![vec![]; hm_size];
                            if use_blockmers {
                                for blockmer in blockmer_info.into_iter() {
                                    let hash = blockmer.kmer % threads as u64;
                                    blockmer_vec_and_canon[hash as usize].push(blockmer);
                                }
                            }

                            // Send both to channels (need to update txs to handle pairs)
                            for (i, (snp_vec, block_vec)) in snpmer_vec_and_canon.into_iter().zip(blockmer_vec_and_canon.into_iter()).enumerate(){
                                txs[i].send((snp_vec, block_vec)).unwrap();
                            }

                            {
                                let mut counter = clone_counter.lock().unwrap();
                                *counter += 1;
                                if *counter % 100000 == 0{
                                    log::info!("Processed {} reads.", counter);
                                }
                            }
                        }
                        Err(_) => {
                            break;
                        }
                    }
                }
                for tx in txs{
                    drop(tx);
                }
            });
        }

        //C: Update bloom filters for both SNPmers and blockmers
        let mut handles = Vec::new();
        for rx in rxs.into_iter(){
            handles.push(thread::spawn(move || {
                // Bloom filters for SNPmers
                let mut snp_filter_canonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut snp_filter_noncanonical = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&42).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut snpmer_map: FxHashMap<u64,[u32;2]> = FxHashMap::default();

                // Bloom filters for blockmers
                let mut block_filter_forward = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&43).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut block_filter_reverse = BloomFilter::with_num_bits((bf_size * 4. * 1_000_000_000. / threads as f64) as usize).seed(&44).expected_items((bf_size * 4. * 1_000_000_000. / 10. / threads as f64) as usize);
                let mut blockmer_map: FxHashMap<u64,[u32;2]> = FxHashMap::default();

                loop{
                    match rx.recv() {
                        Ok(msg) => {
                            let (snpmer_vecs, blockmer_vecs) = msg;

                            // Process SNPmers
                            for kmer_i_canon in snpmer_vecs{
                                let canonical = kmer_i_canon >> 63;
                                let kmer = kmer_i_canon & mask;
                                let kmer_canon = kmer | (1 << 63);
                                if canonical == 1{
                                    let already_present_canon = snp_filter_canonical.insert(&kmer_canon);
                                    let already_present_noncanon = snp_filter_noncanonical.contains(&kmer);
                                    if aggressive_bloom{
                                        if already_present_noncanon && already_present_canon{
                                            snpmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_noncanon{
                                            snpmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                                else{
                                    let already_present_noncanon = snp_filter_noncanonical.insert(&kmer);
                                    let already_present_canon = snp_filter_canonical.contains(&kmer_canon);
                                    if aggressive_bloom{
                                        if already_present_noncanon && already_present_canon{
                                            snpmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_canon{
                                            snpmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                            }

                            // Process blockmers
                            for blockmer in blockmer_vecs{
                                let kmer = blockmer.kmer;

                                if blockmer.is_forward {
                                    let already_present_forward = block_filter_forward.insert(&kmer);
                                    let already_present_reverse = block_filter_reverse.contains(&kmer);
                                    if aggressive_bloom{
                                        if already_present_reverse && already_present_forward{
                                            blockmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_reverse{
                                            blockmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                                else{
                                    let already_present_reverse = block_filter_reverse.insert(&kmer);
                                    let already_present_forward = block_filter_forward.contains(&kmer);
                                    if aggressive_bloom{
                                        if already_present_reverse && already_present_forward{
                                            blockmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                    else{
                                        if already_present_forward{
                                            blockmer_map.insert(kmer, [0,0]);
                                        }
                                    }
                                }
                            }
                        }
                        Err(_) => {
                            log::trace!("Thread finished.");
                            break;
                        }
                    }
                }
                snpmer_map.shrink_to_fit();
                blockmer_map.shrink_to_fit();
                (snpmer_map, blockmer_map)
            }));
        }

        for (map_ind, handle) in handles.into_iter().enumerate() {
            let (snpmer_map, blockmer_map) = handle.join().unwrap();
            bf_snpmer_vec_maps[map_ind] = snpmer_map;
            bf_blockmer_vec_maps[map_ind] = blockmer_map;
        };

    }

    return (bf_snpmer_vec_maps, bf_blockmer_vec_maps);
}

fn second_iteration(
    k: usize,
    l: usize,
    threads: usize,
    args: &Cli,
    bf_snpmer_maps: Vec<FxHashMap<u64, [u32;2]>>,
    bf_blockmer_maps: Vec<FxHashMap<u64, [u32;2]>>) -> (Vec<FxHashMap<u64, [u32;2]>>, Vec<FxHashMap<u64, [u32;2]>>){

    let bf_size = args.bloom_filter_size;
    let mask = !(1 << 63);
    let mut snpmer_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); threads];
    let mut blockmer_vec_maps : Vec<FxHashMap<u64, [u32;2]>> = vec![FxHashMap::default(); threads];
    let minimum_bq = args.minimum_base_quality;

    let num_b = threads/10 + 1;
    let counter = Arc::new(Mutex::new(0));
    let mut rxs = vec![];
    let mut txs_vecs = vec![vec![]; num_b];
    for _ in 0..threads {
        //let (tx, rx) = unbounded();
        let (tx, rx) = bounded(500);
        for i in 1..num_b{
            txs_vecs[i].push(tx.clone());
        }
        txs_vecs[0].push(tx);
        rxs.push(rx);
    }

    //let (tx_head, rx_head1) = unbounded();
    let (tx_head, rx_head1) = bounded(500);
    let mut rx_heads = vec![];
    for _ in 1..num_b{
        let rx_head2 = rx_head1.clone();
        rx_heads.push(rx_head2);
    }
    rx_heads.push(rx_head1);

    assert!(txs_vecs.len() == rx_heads.len());

    let fq_files = args.input_files.clone();
    thread::spawn(move || {
        for fq_file in fq_files{
            let bufreader = BufReader::new(std::fs::File::open(fq_file).expect("valid path"));
            let mut reader = needletail::parse_fastx_reader(bufreader).expect("valid path");
            while let Some(record) = reader.next() {
                let rec = record.expect("Error reading record");
                    let rec_id_string = std::str::from_utf8(rec.id()).unwrap();
                    let last_field = rec_id_string.split_whitespace().last().unwrap_or(""); 
                    let seq;
                    let qualities; 
                    if last_field == "rc"{
                        seq = reverse_complement(&rec.seq().to_vec());
                        qualities = rec.qual().map(|q| q.iter().rev().cloned().collect());
                    }
                    else{
                        seq = rec.seq().to_vec();
                        qualities = rec.qual().map(Vec::from);
                    }
                tx_head.send((seq,qualities)).unwrap();
            }
        }
        drop(tx_head);
        log::debug!("Finished reading all reads.");
    });

    //B: Process kmers and send to hash maps
    let use_blockmers = args.use_blockmers;
    for (rx_head, txs) in rx_heads.into_iter().zip(txs_vecs.into_iter()){
        let clone_counter = Arc::clone(&counter);
        thread::spawn(move || {
            loop{
                match rx_head.recv() {
                    Ok((seq,qualities)) => {
                        // Extract both SNPmers and blockmers
                        let split_kmer_info = seeding::split_kmer_mid(seq.clone(), qualities.clone(), k, minimum_bq);

                        // Hash and distribute SNPmers
                        let mut snpmer_vec_and_canon = vec![vec![]; threads];
                        for kmer_i_and_canon in split_kmer_info.into_iter() {
                            let kmer = kmer_i_and_canon & mask;
                            let hash = kmer % threads as u64;
                            snpmer_vec_and_canon[hash as usize].push(kmer_i_and_canon);
                        }

                        let mut blockmer_vec_and_canon = vec![vec![]; threads];
                        if use_blockmers {
                            let blockmer_info = seeding::blockmer_kmers(seq, qualities, k, l, minimum_bq);
                            // Hash and distribute blockmers
                            for blockmer in blockmer_info.into_iter() {
                                let hash = blockmer.kmer % threads as u64;
                                blockmer_vec_and_canon[hash as usize].push(blockmer);
                            }
                        }

                        // Send both to channels
                        for (i, (snp_vec, block_vec)) in snpmer_vec_and_canon.into_iter().zip(blockmer_vec_and_canon.into_iter()).enumerate(){
                            txs[i].send((snp_vec, block_vec)).unwrap();
                        }

                        {
                            let mut counter = clone_counter.lock().unwrap();
                            *counter += 1;
                            if *counter % 100000 == 0{
                                log::debug!("Processed {} reads.", counter);
                            }
                        }
                    }
                    Err(_) => {
                        break;
                    }
                }
            }
            for tx in txs{
                drop(tx);
            }
        });
    }

    let mut handles = Vec::new();
    for ((rx, my_snpmer_map), my_blockmer_map) in rxs.into_iter().zip(bf_snpmer_maps.into_iter()).zip(bf_blockmer_maps.into_iter()){
        handles.push(thread::spawn(move || {
            let mut my_snpmer_map = my_snpmer_map;
            let mut my_blockmer_map = my_blockmer_map;
            loop{
                match rx.recv() {
                    Ok(msg) => {
                        let (snpmer_vecs, blockmer_vecs) = msg;

                        // Process SNPmers
                        if bf_size > 0.{
                            for kmer_and_canon in snpmer_vecs.into_iter(){
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                if let Some(val) = my_snpmer_map.get_mut(&kmer){
                                    val[canon as usize] += 1;
                                }
                            }
                        }
                        else{
                            for kmer_and_canon in snpmer_vecs.into_iter(){
                                let kmer = kmer_and_canon & mask;
                                let canon = kmer_and_canon >> 63;
                                let val = my_snpmer_map.entry(kmer).or_insert([0,0]);
                                val[canon as usize] += 1
                            }
                        }

                        // Process blockmers
                        if bf_size > 0.{
                            for blockmer in blockmer_vecs.into_iter(){
                                if let Some(val) = my_blockmer_map.get_mut(&blockmer.kmer){
                                    val[blockmer.is_forward as usize] += 1;
                                }
                            }
                        }
                        else{
                            for blockmer in blockmer_vecs.into_iter(){
                                let val = my_blockmer_map.entry(blockmer.kmer).or_insert([0,0]);
                                val[blockmer.is_forward as usize] += 1
                            }
                        }
                    }
                    Err(_) => {
                        log::trace!("Thread finished.");
                        break;
                    }
                }
            }
            my_snpmer_map.shrink_to_fit();
            my_blockmer_map.shrink_to_fit();
            (my_snpmer_map, my_blockmer_map)
        }));
    }

    for (i,handle) in handles.into_iter().enumerate() {
        let (snpmer_map, blockmer_map) = handle.join().unwrap();
        snpmer_vec_maps[i] = snpmer_map;
        blockmer_vec_maps[i] = blockmer_map;
    };

    return (snpmer_vec_maps, blockmer_vec_maps);
}

// Intuitively I think this helps, not sure if we want to use it TODO
pub fn quality_pool(qualities: Vec<u8>) -> Vec<u8>{
    let pool_width = 5;
    let mut pool = Vec::new();
    for i in 0..qualities.len(){
        if i > pool_width/2 && i < qualities.len() - pool_width/2{
            let mut min = 255;
            for j in i-pool_width/2..i+pool_width/2{
                if qualities[j] < min{
                    min = qualities[j];
                }
            }
            pool.push(min);
        }
        else{
            pool.push(qualities[i]);
        }
    }
    return pool;
}