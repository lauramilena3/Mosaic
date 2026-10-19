rule downloadCytoscape:
	output:
		library=os.path.join(workflow.basedir, "tools/cytoscape-3.33.1/cytoscape.min.js"),
		license=os.path.join(workflow.basedir, "tools/cytoscape-3.33.1/LICENSE"),
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	message:
		"Downloading pinned Cytoscape.js for the offline isolate network"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadCytoscape/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p $(dirname {output.library:q})
		wget -q --tries=3 -O {output.library:q} https://cdn.jsdelivr.net/npm/cytoscape@3.33.1/dist/cytoscape.min.js
		printf 'f55947f3daa3bae53209d4b885c195c157f595c225e508a6b382598d9452d6e2  %s\n' {output.library:q} | sha256sum -c -
		wget -q --tries=3 -O {output.license:q} https://cdn.jsdelivr.net/npm/cytoscape@3.33.1/LICENSE
		"""


rule get_SRAToolkit:
	output:
		SRAToolkit_dir=directory("tools/sratoolkit.2.10.0-ubuntu64"),
	message:
		"Downloading SRA toolkit"
	params:
		tools="tools",
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_entrez.yaml",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_SRAToolkit/tot.tsv"
	threads:
		16
	shell:
		"""
		cd {params.tools}
		wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz
		tar -xzf sratoolkit.2.10.0-ubuntu64.tar.gz
		"""

rule downloadContaminants:
	output:
		contaminant_fasta=dirs_dict["CONTAMINANTS_DIR_DB"] +"/{contaminant}.fasta",
		contaminant_dir=temp(directory(dirs_dict["CONTAMINANTS_DIR_DB"] +"/temp_{contaminant}")),
	message:
		"Downloading contaminant genomes"
	params:
		contaminants_dir=dirs_dict["CONTAMINANTS_DIR_DB"],
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_entrez.yaml",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadContaminants/contaminant={contaminant}.tsv"
	threads:
		16
	shell:
		"""
		mkdir -p {output.contaminant_dir}
		cd {output.contaminant_dir}
		wget $(esearch -db "assembly" -query {wildcards.contaminant} | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq | awk -F"/" '{{print $0"/"$NF"_genomic.fna.gz"}}')
		gunzip -f *gz
		cat *fna >> {output.contaminant_fasta}
		"""
		
rule get_VIBRANT:
	output:
		VIBRANT_dir=directory(os.path.join(workflow.basedir, config['vibrant_dir'])),
	message:
		"Downloading VIBRANT"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_VIBRANT/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/AnantharamanLab/VIBRANT
		chmod -R 744 VIBRANT
		cd VIBRANT/databases
		./VIBRANT_setup.py
		"""

rule get_minced:
	output:
		minced_dir=directory(os.path.join(workflow.basedir, config['minced_dir'])),
	message:
		"Downloading and building MinCED"
	conda:
		dirs_dict["ENVS_DIR"] + "/bacterial.yaml"
	params:
		parent=os.path.dirname(os.path.join(workflow.basedir, config['minced_dir'])),
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_minced/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {params.parent:q}
		if [ ! -f {output.minced_dir:q}/Makefile ]; then
			git clone https://github.com/ctSkennerton/minced/ {output.minced_dir:q}
		fi
		make -B -C {output.minced_dir:q} JC="$CONDA_PREFIX/bin/javac" JAR="$CONDA_PREFIX/bin/jar"
		"""

rule get_mmseqs:
	output:
		mmseqs_dir=directory(os.path.join(workflow.basedir, config['mmseqs_dir'])),
		refseq=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna")),
		refseq_taxid=(os.path.join(workflow.basedir,"db/ncbi-taxdump/RefSeqViral.fna.taxidmapping")),
	message:
		"Downloading MMseqs2"
	params:
		taxdump=(os.path.join(workflow.basedir,"db/ncbi-taxdump/")),
	conda:
		dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_mmseqs/tot.tsv"
	threads: 8
	shell:
		"""
		MM_dir={output.mmseqs_dir}
		echo $MM_dir
		if [ ! -d $MM_dir ]
		then
			mkdir -p tools
			cd tools
			git clone https://github.com/soedinglab/MMseqs2.git
			cd MMseqs2
			mkdir build
			cd build
			cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
			make -j {threads}
			make install
		fi
		#download taxdump
		cd ../../../db
		mkdir -p ncbi-taxdump
		cd ncbi-taxdump
		wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
		tar -xzvf taxdump.tar.gz
		#download RefSeqViral
		wget ftp://ftp.ncbi.nlm.nih.gov/blast//db/ref_viruses_rep_genomes.tar.gz
		tar xvzf ref_viruses_rep_genomes.tar.gz
		blastdbcmd -db ref_viruses_rep_genomes -entry all > {output.refseq}
		blastdbcmd -db ref_viruses_rep_genomes -entry all -outfmt "%a %T" > {output.refseq_taxid}
		{output.mmseqs_dir}/build/bin/mmseqs createdb {output.refseq} RefSeqViral.fnaDB
		{output.mmseqs_dir}/build/bin/mmseqs createtaxdb RefSeqViral.fnaDB tmp --ncbi-tax-dump {params.taxdump} --tax-mapping-file {output.refseq_taxid}
		"""

rule get_ALE:
	output:
		ALE_dir=directory(config['ALE_dir']),
	message:
		"Downloading ALE"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_ALE/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/sc932/ALE.git
		cd ALE/src
		make
		"""

rule get_weeSAM:
	output:
		weesam_dir=directory(config['weesam_dir']),
	message:
		"Downloading weesam"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_weeSAM/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/centre-for-virus-research/weeSAM
		"""

rule get_VIGA:
	output:
		VIGA_dir=directory(os.path.join(workflow.basedir, config['viga_dir'])),
		piler_dir=directory(os.path.join(workflow.basedir, config['piler_dir'])),
		# trf_dir=directory(os.path.join(workflow.basedir, config['trf_dir'])),
	message:
		"Downloading MMseqs2"
	# conda:
	# 	dirs_dict["ENVS_DIR"] + "/viga.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_VIGA/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone --depth 1 https://github.com/lauramilena3/viga.git
		chmod 744 viga/create_dbs.sh viga/VIGA.py
		cd viga
		./create_dbs.sh
		cd ..
		wget https://www.drive5.com/pilercr/pilercr1.06.tar.gz --no-check-certificate
		tar -xzvf pilercr1.06.tar.gz
		cd pilercr1.06
		make
		cd ..
		cd TRF
		# wget wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
		# wget http://tandem.bu.edu/irf/downloads/irf307.linux.exe
		# mv trf409.linux64 trf
		# mv irf307.linux.exe irf
		chmod 744 trf irf
		"""

rule downloadCenoteDB:
	output:
		cenote_db=directory(config["cenote_db"]),
	message:
		"Downloading the core Cenote-Taker3 databases"
	conda:
		dirs_dict["ENVS_DIR"] + "/cenote.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadCenoteDB/tot.tsv"
	log:
		dirs_dict["BENCHMARKS"] + "/downloadCenoteDB/tot.log"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		"""
		mkdir -p {output.cenote_db:q}
		get_ct3_dbs -o {output.cenote_db:q} --hmm T --hallmark_tax T --refseq_tax T --mmseqs_cdd T --domain_list T > {log:q} 2>&1
		# The downloader does not propagate all failed download/build commands.
		for database_file in hmmscan_DBs/v3.1.1/Virion_HMMs.h3m hmmscan_DBs/v3.1.1/DNA_rep_HMMs.h3m \\
			hmmscan_DBs/v3.1.1/RDRP_HMMs.h3m hmmscan_DBs/v3.1.1/Useful_Annotation_HMMs.h3m \\
			hmmscan_DBs/v3.1.1/phrogs_for_ct.h3m mmseqs_DBs/ct3_hallmark.taxDB \\
			mmseqs_DBs/refseq_virus_prot_taxDB mmseqs_DBs/CDD viral_cdds_and_pfams_191028.txt; do
			test -s {output.cenote_db:q}/"$database_file"
		done
		date -u '+%Y-%m-%d' > {output.cenote_db:q}/download_date.txt
		"""

rule downloadVirSorterDB:
	output:
		virSorter_dir=directory(config['virSorter_db']),
	message:
		"Downloading VirSorter database"
	threads: 8
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadVirSorterDB/tot.tsv"
	params:
		virSorter_db="db/VirSorter"
	shell:
		"""
		#git clone https://github.com/jiarong/VirSorter2.git {output.virSorter_dir}
		#cd {output.virSorter_dir}
		#pip install .
		virsorter setup -d {output.virSorter_dir} -j {threads}
		#mkdir {output.virSorter_dir}
		"""

rule downloadIphopDB:
	output:
		iphop_db=directory(config['iphop_db']),
	message:
		"Downloading iphop database"
	threads: 1
	params:
		db_dir="db/iphop_db/"
	conda:
		dirs_dict["ENVS_DIR"] + "/env2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadIphopDB/tot.tsv"
	shell:
		"""
		mkdir -p {params.db_dir}
		iphop download --db_dir {params.db_dir} --db_version iPHoP_db_Jun25_rw --split --no_prompt
		iphop download --db_dir {output.iphop_db} --full_verify
		"""

# Serial component setup avoids concurrent writes to the database-local CONFIG.
DRAM_COMPONENTS=["kofam_hmm", "kofam_ko_list", "uniref", "pfam", "pfam_hmm", "dbcan",
	"dbcan_fam_activities", "dbcan_subfam_ec", "vogdb", "vog_annotations", "viral", "peptidase",
	"genome_summary_form", "module_step_form", "function_heatmap_form", "amg_database", "etc_module_database"]
DRAM_SETUP_DIR=config['DRAM_db'].rstrip("/") + "/.mosaic_setup"

def input_previous_dram_component(wildcards):
	index=DRAM_COMPONENTS.index(wildcards.dram_component)
	return [DRAM_SETUP_DIR + "/" + DRAM_COMPONENTS[index - 1] + ".json"] if index else []

rule prepareDRAMcomponent:
	input:
		input_previous_dram_component,
	output:
		checkpoint=DRAM_SETUP_DIR + "/{dram_component}.json",
	params:
		db_dir=config['DRAM_db'].rstrip("/"),
	wildcard_constraints:
		dram_component="|".join(DRAM_COMPONENTS),
	message:
		"Preparing DRAM database component {wildcards.dram_component}"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/prepareDRAMcomponent/component={dram_component}.tsv"
	log:
		config['DRAM_db'].rstrip("/") + ".{dram_component}.log"
	shell:
		r"""
		mkdir -p "$(dirname {log:q})"
		(
		python - {params.db_dir:q} {wildcards.dram_component:q} {output.checkpoint:q} {threads} <<-'PYTHON'
		import ast
		from datetime import datetime
		import gzip
		import importlib
		import importlib.util
		import json
		import logging
		import os
		from pathlib import Path
		import re
		import shutil
		import sys
		import tarfile
		from urllib.request import Request, urlopen

		db_dir=Path(sys.argv[1]).resolve()
		component=sys.argv[2]
		checkpoint=Path(sys.argv[3])
		state_dir=db_dir / ".mosaic_setup"
		state_dir.mkdir(parents=True, exist_ok=True)
		configuration=state_dir / "CONFIG"
		package=Path(importlib.util.find_spec("mag_annotator").origin).parent
		if configuration.exists():
		    config=json.loads(configuration.read_text())
		else:
		    config=json.loads((package / "CONFIG").read_text())
		    # Adopt only existing files from this database, not another installation.
		    for section in ["search_databases", "database_descriptions", "dram_sheets"]:
		        for name, value in config[section].items():
		            candidate=Path(value).resolve() if value else None
		            if not candidate or db_dir not in candidate.parents or not candidate.is_file() or not candidate.stat().st_size:
		                config[section][name]=None
		    description=Path(config["description_db"]).resolve() if config.get("description_db") else None
		    if not description or db_dir not in description.parents or not description.is_file():
		        config["description_db"]=None
		    configuration.write_text(json.dumps(config, indent=2))
		section=next((name for name in ["search_databases", "database_descriptions", "dram_sheets"]
		    if component in config[name]), None)
		if section is None:
		    raise RuntimeError("Unknown DRAM configuration component: " + component)
		value=config[section].get(component)
		def component_artifacts(value):
		    path=Path(value).resolve()
		    files=[path]
		    if component in ["kofam_hmm", "dbcan", "vogdb"]:
		        files += [Path(str(path) + suffix) for suffix in [".h3f", ".h3i", ".h3m", ".h3p"]]
		    if component in ["uniref", "pfam", "viral", "peptidase"]:
		        files += [Path(str(path) + suffix) for suffix in [".dbtype", ".index", "_h", ".idx", ".idx.dbtype", ".idx.index"]]
		    return files
		def valid_artifacts(files):
		    if not files or not all(db_dir in path.parents and path.is_file() and path.stat().st_size for path in files):
		        return False
		    path=files[0]
		    try:
		        opener=gzip.open if path.suffix == ".gz" else open
		        with opener(path, "rb") as handle:
		            prefix=handle.read(8192).lstrip()
		        if component in ["kofam_hmm", "dbcan", "vogdb"]:
		            return prefix.startswith(b"HMMER3/")
		        if component in ["dbcan_fam_activities", "dbcan_subfam_ec"]:
		            rows=[line for line in prefix.splitlines() if line.strip() and not line.startswith(b"#")]
		            return bool(rows) and not prefix.startswith(b"<") and all(b"\t" in line for line in rows[:3])
		        if component == "vog_annotations":
		            return all(field in prefix.splitlines()[0].split(b"\t") for field in [b"FunctionalCategory", b"ConsensusFunctionalDescription"])
		        return bool(prefix) and not prefix.lower().startswith((b"<!doctype html", b"<html"))
		    except (OSError, EOFError, IndexError):
		        return False
		artifacts=component_artifacts(value) if value else []
		if valid_artifacts(artifacts):
		    checkpoint.write_text(json.dumps(dict(component=component, path=value, reused=True,
		        artifacts=[str(path) for path in artifacts]), indent=2))
		    print("Reusing completed DRAM component: " + component, flush=True)
		    sys.exit(0)

		# Patch the active environment, never a server-specific Conda path.
		source=Path(importlib.util.find_spec("mag_annotator.database_processing").origin)
		original=source.read_text()
		base="https://pro.unl.edu/dbCAN2/download/Databases/"
		old_base=r"https?://bcb\.unl\.edu/dbCAN2/download/"
		updated=re.sub(old_base + r"(?:Databases/)?(?:V\{{version\}}/)?(?=dbCAN-HMMdb-V\{{version\}}\.txt)",
		    base + "V{{version}}/", original)
		updated=re.sub(old_base + "Databases/", base, updated)
		updated=updated.replace("http://fileshare.csb.univie.ac.at/vog/", "https://fileshare.csb.univie.ac.at/vog/")
		old_vog="merge_files(glob(path.join(hmm_dir, 'VOG*.hmm')), vog_hmms)"
		new_vog="vog_files=sorted(glob(path.join(hmm_dir, '**', 'VOG*.hmm'), recursive=True))\n    if not vog_files:\n        raise ValueError('No VOG HMM profiles found in the downloaded archive')\n    merge_files(vog_files, vog_hmms)"
		updated=updated.replace(old_vog, new_vog)
		functions=["download_dbcan", "download_dbcan_fam_activities", "download_dbcan_subfam_ec"]
		tree=ast.parse(updated, filename=str(source))
		available=[node.name for node in tree.body if isinstance(node, ast.FunctionDef)]
		if re.search(old_base, updated) or base not in updated or not all(name in available for name in functions) or new_vog not in updated:
		    raise RuntimeError("Unexpected DRAM dbCAN/VOG download code; no patch was written. Review the installed DRAM version.")
		if updated != original:
		    backup=source.with_name(source.name + ".before_mosaic_database_fixes")
		    if not backup.exists():
		        shutil.copy2(source, backup)
		    source.write_text(updated)
		    print("Patched DRAM dbCAN URLs and nested VOG archive handling in " + str(source), flush=True)
		else:
		    print("DRAM dbCAN/VOG fixes are already patched.", flush=True)

		# Check the selected dbCAN/VOG input before starting its real download.
		database_processing=importlib.import_module("mag_annotator.database_processing")
		def check_download(url, output_file, logger, verbose=True):
		    with urlopen(Request(url, headers={{"Range": "bytes=0-8191"}}), timeout=60) as response:
		        if component == "vogdb":
		            with tarfile.open(fileobj=response, mode="r|gz") as archive:
		                for member in archive:
		                    if member.isfile() and Path(member.name).match("VOG*.hmm"):
		                        if not archive.extractfile(member).read(16).startswith(b"HMMER3/"):
		                            raise RuntimeError("Invalid VOG HMM profile in " + url)
		                        print("Checked nested VOG archive: " + member.name, flush=True)
		                        return
		                raise RuntimeError("No VOG HMM profiles found in " + url)
		        if component == "vog_annotations":
		            with gzip.GzipFile(fileobj=response) as handle:
		                prefix=handle.read(4096).decode("utf-8")
		            if not all(field in prefix.splitlines()[0].split("\t") for field in ["FunctionalCategory", "ConsensusFunctionalDescription"]):
		                raise RuntimeError("Unexpected VOG annotation columns: " + url)
		            print("Checked VOG annotations: " + url, flush=True)
		            return
		        prefix=response.read(8192).decode("utf-8", errors="replace").lstrip()
		    if Path(output_file).name.startswith("dbCAN-HMMdb-"):
		        valid=prefix.startswith("HMMER3/")
		    else:
		        rows=[line for line in prefix.splitlines() if line.strip() and not line.startswith("#")]
		        valid=bool(rows) and not prefix.startswith("<") and all("\t" in line for line in rows[:3])
		    if not valid:
		        raise RuntimeError("DRAM dbCAN download is not the expected HMM/table (possibly HTML): " + url)
		    print("Checked dbCAN download: " + url, flush=True)
		download=getattr(database_processing, "download_" + component)
		if component in ["dbcan", "dbcan_fam_activities", "dbcan_subfam_ec", "vogdb", "vog_annotations"]:
		    real_download=database_processing.download_file
		    database_processing.download_file=check_download
		    try:
		        download(".")
		    finally:
		        database_processing.download_file=real_download

		component_dir=db_dir / "components" / component
		saved_configuration=configuration.read_bytes()
		if component_dir.exists():
		    if not (component_dir / ".mosaic_owned").is_file():
		        raise RuntimeError("Refusing to replace a component directory not created by MOSAIC: " + str(component_dir))
		    shutil.rmtree(component_dir)
		component_dir.mkdir(parents=True)
		(component_dir / ".mosaic_owned").touch()
		raw_dir=component_dir / "raw"
		raw_dir.mkdir()
		logger=logging.getLogger("mosaic_dram")
		logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(message)s")
		try:
		    raw=download(str(raw_dir), logger)
		    process=getattr(database_processing, "process_" + component, None)
		    if process:
		        paths=process(raw, str(component_dir), logger, threads=int(sys.argv[4]))
		    else:
		        final=component_dir / Path(raw).name
		        shutil.move(raw, final)
		        paths={{component: str(final)}}
		    artifacts=component_artifacts(paths[component])
		    if not valid_artifacts(artifacts):
		        raise RuntimeError("DRAM component did not produce valid data/index files: " + component)
		    if config.get("setup_info") is None:
		        config["setup_info"]={{}}
		    for name, value in paths.items():
		        destination=next(key for key in ["search_databases", "database_descriptions", "dram_sheets"] if name in config[key])
		        config[destination][name]=str(Path(value).resolve())
		        config.setdefault("setup_info", {{}})[name]=dict(name=name,
		            Origin="Downloaded and processed by MOSAIC using DRAM", **{{"Download time": datetime.now().isoformat()}})
		    # Only paths from a completely successful component enter the shared CONFIG.
		    config["description_db"]=None
		    pending=configuration.with_suffix(".pending")
		    pending.write_text(json.dumps(config, indent=2))
		    os.replace(pending, configuration)
		    for name in ["raw", "tmp", "kofam_profiles", "vogdb_hmms"]:
		        scratch=component_dir / name
		        if scratch.is_dir():
		            shutil.rmtree(scratch)
		    checkpoint.write_text(json.dumps(dict(component=component, paths=paths, reused=False), indent=2))
		except BaseException:
		    # Never remove the database root or another completed component.
		    shutil.rmtree(component_dir)
		    pending=configuration.with_suffix(".pending")
		    pending.write_bytes(saved_configuration)
		    os.replace(pending, configuration)
		    raise
		PYTHON
		) 2>&1 | tee {log:q}
		"""

rule downloadDRAMDB:
	input:
		components=[DRAM_SETUP_DIR + "/" + name + ".json" for name in DRAM_COMPONENTS],
	output:
		ready=DRAM_SETUP_DIR + "/complete.json",
	params:
		db_dir=config['DRAM_db'].rstrip("/"),
	conda:
		dirs_dict["ENVS_DIR"] + "/vir2.yaml"
	message:
		"Finalizing the DRAM description database and readiness record"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadDRAMDB/tot.tsv"
	log:
		config['DRAM_db'].rstrip("/") + ".log"
	threads: 4
	shell:
		r"""
		mkdir -p "$(dirname {log:q})"
		python - {params.db_dir:q} {output.ready:q} <<-'PYTHON' 2>&1 | tee {log:q}
		import json
		import logging
		import os
		from pathlib import Path
		import sqlite3
		import sys
		from mag_annotator.database_handler import DatabaseHandler

		db_dir=Path(sys.argv[1]).resolve()
		configuration=db_dir / ".mosaic_setup" / "CONFIG"
		config=json.loads(configuration.read_text())
		logger=logging.getLogger("mosaic_dram")
		logging.basicConfig(level=logging.INFO)
		previous=Path(config["description_db"]) if config.get("description_db") else None
		valid=False
		if previous and previous.is_file():
		    try:
		        with sqlite3.connect(previous.as_uri() + "?mode=ro", uri=True) as database:
		            tables={{row[0] for row in database.execute("SELECT name FROM sqlite_master WHERE type='table'")}}
		            valid=database.execute("PRAGMA quick_check").fetchone()[0] == "ok" and all(
		                name + "_description" in tables and database.execute(
		                    "SELECT 1 FROM " + name + "_description LIMIT 1").fetchone() is not None
		                for name in ["uniref", "pfam", "dbcan", "vogdb", "viral", "peptidase"])
		    except sqlite3.DatabaseError:
		        valid=False
		if not valid:
		    pending=db_dir / ".mosaic_setup" / "description_db.sqlite.pending"
		    if pending.exists():
		        pending.unlink()
		    handler=DatabaseHandler(logger, config_loc=str(configuration))
		    try:
		        handler.populate_description_db(str(pending), update_config=False)
		        handler.session.close()
		        handler.session.bind.dispose()
		        final=db_dir / "description_db.sqlite"
		        os.replace(pending, final)
		        handler.config["description_db"]=str(final)
		        config=handler.config
		    except BaseException:
		        if pending.exists():
		            pending.unlink()
		        raise
		    temporary=configuration.with_suffix(".pending")
		    temporary.write_text(json.dumps(config, indent=2))
		    os.replace(temporary, configuration)
		Path(sys.argv[2]).write_text(json.dumps(dict(configuration=str(configuration), description_db=config["description_db"]), indent=2))
		print("DRAM database ready: " + str(db_dir))
		PYTHON
		"""

rule downloadCheckvDB:
	output:
		checkv_db=directory(config['checkv_db']),
	message:
		"Downloading CheckV database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadCheckvDB/tot.tsv"
	shell:
		"""
		checkv download_database ./db
		"""

rule downloadCheckMDB:
	output:
		checkm_db=directory(config['checkm_db']),
		checkm_tar=temp("checkm_data_2015_01_16.tar.gz")
	message:
		"Downloading CheckM database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadCheckMDB/tot.tsv"
	shell:
		"""
		wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
		mkdir -p {output.checkm_db}
		tar xzvf {output.checkm_tar} -C {output.checkm_db}
		checkm data setRoot {output.checkm_db}
		"""

rule downloadGtdbtk_db:
	output:
		gtdbtk_db=directory(config['gtdbtk_db']),
	message:
		"Downloading GTDB-Tk database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadGtdbtk_db/tot.tsv"
	shell:
		"""
		wget https://data.gtdb.ecogenomic.org/releases/release214/214.0/auxillary_files/gtdbtk_r214_data.tar.gz
		mkdir -p {output.gtdbtk_db}
		tar xzvf gtdbtk_r214_data.tar.gz -C {output.gtdbtk_db}
		"""
		

rule getKrakenTools:
	output:
		kraken_tools=directory(config['kraken_tools']),
	message:
		"Downloading KrakenTools"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getKrakenTools/tot.tsv"
	shell:
		"""
		mkdir -p tools
		cd tools
		git clone https://github.com/jenniferlu717/KrakenTools
		chmod 777 KrakenTools/*
		"""

rule downloadPhaBOXDB:
	output:
		phabox_db=directory(config["phabox_db"]),
		archive=temp(os.path.normpath(config["phabox_db"]) + ".zip"),
	params:
		url=config["phabox_db_url"],
		sha256=config.get("phabox_db_sha256", ""),
	message:
		"Downloading and extracting the PhaBOX v2.2 reference database"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadPhaBOXDB/tot.tsv"
	shell:
		"""
		mkdir -p "$(dirname {output.archive:q})"
		wget --continue --tries=3 --timeout=60 --output-document={output.archive:q} {params.url:q}
		python - {output.archive:q} {output.phabox_db:q} {params.url:q} {params.sha256:q} <<-'PYTHON'
		import csv
		from datetime import datetime, timezone
		import hashlib
		from pathlib import Path, PurePosixPath
		import shutil
		import sys
		import tempfile
		import zipfile
		archive_path, target, source_url, expected_sha256 = sys.argv[1:]
		target = Path(target)
		checksum = hashlib.sha256()
		with open(archive_path, "rb") as handle:
		    for block in iter(lambda: handle.read(1024 * 1024), b""):
		        checksum.update(block)
		if expected_sha256 and checksum.hexdigest() != expected_sha256:
		    raise ValueError("PhaBOX database checksum does not match the configured release")
		with tempfile.TemporaryDirectory(prefix="phabox_download_", dir=target.parent) as temporary:
		    staging = Path(temporary)
		    with zipfile.ZipFile(archive_path) as archive:
		        for member in archive.infolist():
		            path = PurePosixPath(member.filename)
		            if path.is_absolute() or ".." in path.parts:
		                raise ValueError("Unsafe path in the PhaBOX database archive")
		        archive.extractall(staging)
		    database = next(staging.rglob("RefVirus.dmnd")).parent
		    for filename in ["RefVirus.faa", "RefVirus.csv", "database_aai.tsv", "nodes.csv", "names.csv", "taxid.csv", "ProkaryoticGroup.pkl", "RefVirus_anno.pkl"]:
		        if not (database / filename).is_file():
		            raise FileNotFoundError("Missing PhaBOX taxonomy database file: " + filename)
		    target.mkdir(parents=True, exist_ok=True)
		    for path in database.iterdir():
		        shutil.move(str(path), target / path.name)
		with open(target / "download_info.tsv", "w", newline="") as handle:
		    writer = csv.writer(handle, delimiter="\\t")
		    writer.writerow(["database_version", "downloaded_utc", "source_url", "sha256"])
		    writer.writerow(["2.2", datetime.now(timezone.utc).isoformat(), source_url, checksum.hexdigest()])
		PYTHON
		"""

rule downloadKrakenDB:
	output:
		kraken_db=directory(config['kraken_db']),
		kraken_tar=temp("k2_pluspfp_08_GB_20260626.tar.gz")
	message:
		"Downloading Kraken database"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadKrakenDB/tot.tsv"
	shell:
		"""
		wget https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08_GB_20260626.tar.gz
		mkdir -p {output.kraken_db:q}
		tar -xvf {output.kraken_tar} -C {output.kraken_db}
		"""

rule downloadSourmashRocksDB:
	output:
		rocksdb=directory(config["sourmash_rocksdb"]),
		archive=temp(config["sourmash_rocksdb"] + ".tar.gz"),
	params:
		url=config["sourmash_rocksdb_url"],
	message:
		"Downloading the GTDB RS226 k31 Sourmash RocksDB index"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadSourmashRocksDB/tot.tsv"
	threads: 1
	resources:
		mem_mb=1000,
	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname {output.archive:q})" {output.rocksdb:q}
		wget --continue --tries=3 --output-document={output.archive:q} {params.url:q}
		tar -xzf {output.archive:q} --strip-components=1 -C {output.rocksdb:q}
		"""

rule downloadSourmashTaxonomy:
	output:
		taxonomy=config["sourmash_tax"],
	params:
		url=config["sourmash_tax_url"],
	message:
		"Downloading the matching GTDB RS226 Sourmash taxonomy"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadSourmashTaxonomy/tot.tsv"
	threads: 1
	resources:
		mem_mb=1000,
	shell:
		"""
		set -euo pipefail
		mkdir -p "$(dirname {output.taxonomy:q})"
		wget --continue --tries=3 --output-document={output.taxonomy:q} {params.url:q}
		"""

rule downloadKrakenUniqDB:
	output:
		krakenuniq_db=directory(config['krakenUniq_db']),
		krakenuniq_tar=temp("kuniq_standard_minus_kdb.20220616.tgz")
	message:
		"Downloading Kraken database"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"] + "/env4.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadKrakenUniqDB/tot.tsv"
	shell:
		"""
		mkdir {output.krakenuniq_db}
		cd {output.krakenuniq_db}
		wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/kuniq_standard_minus_kdb.20220616.tgz
		wget https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2022-06-16-STANDARD/database.kdb 
		tar -xvf {output.krakenuniq_tar} -C {output.krakenuniq_db}
		"""


rule installBracken:
	output:
		bracken_dir=directory(config['bracken_dir']),
	message:
		"Downloading Braken"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/installBracken/tot.tsv"
	shell:
		"""
		mkdir -p tools
		cd tools
    	git clone https://github.com/jenniferlu717/Bracken
		cd Bracken
		bash install_bracken.sh
		"""

rule buildBrackenDB:
	input:
		bracken_dir=(config['bracken_dir']),
		kraken_db=config['kraken_db'],
	output:
		bracken_checkpoint=config['kraken_db'] + "../bracken_db_ckeckpoint.txt",
	message:
		"Building Braken database"
	threads: 144
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBrackenDB/tot.tsv"
	shell:
		"""
    	bracken-build -d {input.kraken_db} -t {threads} -k 35 -l 150
		echo {input.kraken_db} > {output.bracken_checkpoint}
		"""

rule buildBrackenUniqDB:
	input:
		bracken_dir=config['bracken_dir'],
		krakenuniq_db=config['krakenUniq_db'],
	output:
		brackenuniq_checkpoint=config['krakenUniq_db'] + "_brackenuniq_db_ckeckpoint.txt",
	message:
		"Building BrakenUniq database"
	threads: 32
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/buildBrackenUniqDB/tot.tsv"
	shell:
		"""
    	{input.bracken_dir}/bracken-build -d {input.krakenuniq_db} -t {threads} -k 31 -l 150 -y krakenuniq
		touch {output.brackenuniq_checkpoint}
		"""

rule downloadGenomadDB:
	output:
		genomad_db=directory(config['genomad_db']),
	message:
		"Downloading geNomad database"
	params:
		db_dir="db/"
	conda:
		dirs_dict["ENVS_DIR"] + "/env6.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadGenomadDB/tot.tsv"
	shell:
		"""
		genomad download-database {params.db_dir}
		"""

rule downloadTaxmyphageDB:
	output:
		taxmyphage_db=directory(config['taxmyphage_db']),
	message:
		"Downloading taxmyphage database"
	params:
		db_dir="db/"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadTaxmyphageDB/tot.tsv"
	shell:
		"""
		taxmyphage install -db {output.taxmyphage_db}
		cd {output.taxmyphage_db} 
		# wget https://ictv.global/sites/default/files/VMR/VMR_MSL39_v1.xlsx
		# mv VMR_MSL39_v1.xlsx VMR.xlsx 
		"""

rule downloadPharokkaDB:
	output:
		pharokka_db=directory(config['pharokka_db']),
	message:
		"Downloading pharokka database"
	params:
		db_dir="db/"
	conda:
		dirs_dict["ENVS_DIR"] + "/env7.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadPharokkaDB/tot.tsv"
	shell:
		"""
		install_databases.py -o {output.pharokka_db}
		"""

rule downloadPhyntenyDB:
	output:
		phynteny_db=directory(config['phynteny_db']),
	message:
		"Downloading phynteny models"
	conda:
		dirs_dict["ENVS_DIR"] + "/phynteny.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadPhyntenyDB/tot.tsv"
	shell:
		"""
		install_models -o {output.phynteny_db}
		"""
		
rule downloadKrakenDB_human:
	output:
		kraken_db_human=directory(config['kraken_db_human']),
	message:
		"Downloading human Kraken database"
	threads: 4
	conda:
		dirs_dict["ENVS_DIR"] + "/env1_kraken.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadKrakenDB_human/tot.tsv"
	shell:
		"""
		kraken2-build --download-library human --db {output.kraken_db_human} --threads {threads} --use-ftp
		kraken2-build --download-taxonomy --db {output.kraken_db_human}
		kraken2-build --build --db {output.kraken_db_human} --threads {threads}
		kraken2-build --clean --db {output.kraken_db_human}
		"""

rule downloadPholdDB:
	output:
		database=directory(config["phold_db"]),
	message:
		"Downloading the Phold database and protein-structure prediction model"
	conda:
		dirs_dict["ENVS_DIR"] + "/phold.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadPholdDB/tot.tsv"
	threads: 8
	resources:
		mem_mb=int(config.get("phold_mem_mb", 32768)),
	shell:
		"""
		phold install -d {output.database:q} -t {threads}
		"""

rule downloadVcontact2Files:
	output:
		archive="db/vcontact2/GenomesDB_Aug_2026.tar.gz",
		gene2genome_millard="db/vcontact2/GenomesDB_Aug_2026_vConTACT2_gene_to_genome.csv",
		vcontact_aa_millard="db/vcontact2/GenomesDB_Aug_2026_vConTACT2_proteins.faa",
	params:
		url=config["millard_genomesdb_url"],
	message:
		"Downloading Millard GenomesDB and preparing vConTACT2 reference inputs"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadVcontact2Files/tot.tsv"
	threads: 1
	shell:
		r"""
		set -euo pipefail
		mkdir -p db/vcontact2
		wget --continue --tries=3 --output-document={output.archive:q} {params.url:q}
		python - {output.archive:q} {output.gene2genome_millard:q} {output.vcontact_aa_millard:q} <<'PYCODE'
import csv
import sys
import tarfile

with tarfile.open(sys.argv[1], "r|gz") as archive, \
        open(sys.argv[2], "w", newline="") as mapping, \
        open(sys.argv[3], "wb") as proteins:
    writer = csv.writer(mapping)
    writer.writerow(["protein_id", "contig_id", "keywords"])
    count = 0
    for member in archive:
        parts = member.name.split("/")
        if not member.isfile() or len(parts) != 3 or parts[0] != "GenomesDB" or not parts[2].endswith(".faa"):
            continue
        with archive.extractfile(member) as source:
            last_line = b""
            for line in source:
                if line.startswith(b">"):
                    protein_id = line[1:].split(None, 1)[0].decode("utf-8")
                    writer.writerow([protein_id, parts[1], "none"])
                    count += 1
                proteins.write(line)
                last_line = line
            if last_line and not last_line.endswith(b"\n"):
                proteins.write(b"\n")
    if count == 0:
        raise ValueError("No protein sequences found in the Millard GenomesDB archive")
PYCODE
		"""

rule downloadMetaVR:
	output:
		archive=os.path.join(config["METAVR_db"], "METAVR_blastdb.tar.zst"),
		blastdb=directory(os.path.join(config["METAVR_db"], "METAVR_UViG_blastdb")),
		metadata=os.path.join(config["METAVR_db"], "METAVR_main_table.parquet"),
		sources=os.path.join(config["METAVR_db"], "IMG_full_metadata.tsv.gz"),
	params:
		url=config["metavr_download_url"].rstrip("/"),
		db_dir=config["METAVR_db"],
		blastdb=os.path.join(config["METAVR_db"], "METAVR_UViG_blastdb", "METAVR_UViG.blastdb"),
	message:
		"Downloading MetaVR BLAST database and metadata from NERSC"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadMetaVR/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		"""
		set -euo pipefail
		mkdir -p {params.db_dir:q}
		wget --continue --tries=3 --output-document={output.archive:q} {params.url:q}/METAVR_blastdb.tar.zst
		wget --continue --tries=3 --output-document={output.metadata:q} {params.url:q}/METAVR_main_table.parquet
		wget --continue --tries=3 --output-document={output.sources:q} {params.url:q}/IMG_full_metadata.tsv.gz
		gzip -t {output.sources:q}
		tar --use-compress-program=zstd -xf {output.archive:q} -C {params.db_dir:q}
		blastdbcmd -db {params.blastdb:q} -info
		"""

rule downloadMetaVRNucleotideFasta:
	output:
		fasta=config["METAVR_reference_fasta"],
	params:
		url=config["metavr_download_url"].rstrip("/") + "/METAVR.fna.bgz",
	message:
		"Downloading the MetaVR nucleotide FASTA for isolate-relative extraction"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadMetaVRNucleotideFasta/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		"""
		set -euo pipefail
		mkdir -p $(dirname {output.fasta:q})
		wget --continue --tries=3 --output-document={output.fasta:q}.bgz {params.url:q}
		gzip -t {output.fasta:q}.bgz
		gzip -dc {output.fasta:q}.bgz > {output.fasta:q}
		test -s {output.fasta:q}
		"""

rule downloadMetaVRProteinFasta:
	output:
		fasta=config["METAVR_protein_db"],
	params:
		url=config["metavr_download_url"].rstrip("/") + "/METAVR.faa.bgz",
	message:
		"Downloading and indexing the MetaVR protein FASTA for isolate searches"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadMetaVRProteinFasta/tot.tsv"
	threads: 1
	resources:
		mem_mb=8000,
	shell:
		"""
		set -euo pipefail
		mkdir -p $(dirname {output.fasta:q})
		wget --continue --tries=3 --output-document={output.fasta:q}.bgz {params.url:q}
		gzip -t {output.fasta:q}.bgz
		gzip -dc {output.fasta:q}.bgz > {output.fasta:q}
		test -s {output.fasta:q}
		makeblastdb -in {output.fasta:q} -dbtype prot
		"""

rule downloadRefSeqViral:
	output:
		fasta="db/RefSeqViral/RefSeq_viral.fasta",
		manifest="db/RefSeqViral/RefSeq_viral.download_info.tsv",
	params:
		url=config.get("refseq_viral_download_url", "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"),
	message:
		"Downloading and dating the viral RefSeq nucleotide database"
	conda:
		dirs_dict["ENVS_DIR"] + "/env5.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadRefSeqViral/tot.tsv"
	threads: 1
	shell:
		"""
		set -euo pipefail
		mkdir -p db/RefSeqViral
		stamp=$(date -u +%Y-%m-%dT%H%M%SZ)
		stage=$(mktemp -d db/RefSeqViral/.download_XXXXXX)
		trap 'rm -r -- "$stage"' EXIT
		wget --tries=3 --server-response --output-document="$stage/source.fna.gz" \
			{params.url:q} 2> "$stage/http_headers.txt"
		gzip -t "$stage/source.fna.gz"
		gzip -dc "$stage/source.fna.gz" > "$stage/RefSeq_viral.fasta"
		test -s "$stage/RefSeq_viral.fasta"
		makeblastdb -in "$stage/RefSeq_viral.fasta" -dbtype nucl -parse_seqids
		test -s "$stage/RefSeq_viral.fasta.nhr"
		test -s "$stage/RefSeq_viral.fasta.nin"
		test -s "$stage/RefSeq_viral.fasta.nsq"
		printf 'download_started_utc\t%s\n' "$stamp" > "$stage/download_info.tsv"
		printf 'download_completed_utc\t%s\n' "$(date -u +%Y-%m-%dT%H%M%SZ)" >> "$stage/download_info.tsv"
		printf 'source_url\t%s\n' {params.url:q} >> "$stage/download_info.tsv"
		printf 'archive_sha256\t%s\n' "$(sha256sum "$stage/source.fna.gz" | cut -d' ' -f1)" >> "$stage/download_info.tsv"
		printf 'fasta_sha256\t%s\n' "$(sha256sum "$stage/RefSeq_viral.fasta" | cut -d' ' -f1)" >> "$stage/download_info.tsv"
		for existing in {output.fasta:q} {output.manifest:q} db/RefSeqViral/RefSeq_viral.fasta.n*; do
			if [ -e "$existing" ] && [ ! -L "$existing" ]; then
				printf 'Refusing to replace existing non-symlink: %s\n' "$existing" >&2
				exit 1
			fi
		done
		snapshot="db/RefSeqViral/$stamp"
		if [ -e "$snapshot" ]; then
			printf 'Dated RefSeq directory already exists: %s\n' "$snapshot" >&2
			exit 1
		fi
		mv "$stage" "$snapshot"
		trap - EXIT
		for index in "$snapshot"/RefSeq_viral.fasta.n*; do
			[ -f "$index" ] || continue
			name=$(basename "$index")
			ln -sfn "$stamp/$name" "db/RefSeqViral/$name"
		done
		ln -sfn "$stamp/RefSeq_viral.fasta" {output.fasta:q}
		ln -sfn "$stamp/download_info.tsv" {output.manifest:q}
		"""

rule downloadBLASTviralProteins:
	output:
		blast=(os.path.join(workflow.basedir,"db/ncbi/NCBI_viral_proteins.faa")),
	message:
		"Downloading RefSeq viral proteins for blast annotation"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_blast_download.yaml",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadBLASTviralProteins/tot.tsv"
	shell:
		"""
		esearch -db "protein" -query "txid10239[Organism:exp] AND (viruses[filter] AND refseq[filter])" \
			| efetch -format fasta > {output.blast}
		makeblastdb -in {output.blast} -dbtype prot
		"""
		
rule getClusterONE:
	output:
		clusterONE_dir=directory(config["clusterONE_dir"]),
	message:
		"Downloading clusterONE"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/getClusterONE/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.clusterONE_dir}
		wget http://www.paccanarolab.org/static_content/clusterone/cluster_one-1.0.jar --no-check-certificate
		mv cluster_one-1.0.jar {output.clusterONE_dir}
		chmod 744 {output.clusterONE_dir}/cluster_one-1.0.jar
		"""
rule downloadCanu:
	output:
		canu_dir=directory(config['canu_dir']),
	message:
		"Installing Canu assembler"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadCanu/tot.tsv"
	threads: 1
	shell:
		"""
		if [ ! -d {output.canu_dir} ]
		then
			if [ {config[operating_system]} == "macOs" ]
			then
				mkdir -p tools
				curl -OL https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Darwin-amd64.tar.xz
			else
				mkdir -p tools
				curl -OL https://github.com/marbl/canu/releases/download/v2.0/canu-2.0.Linux-amd64.tar.xz
			fi
		fi
		tar -xJf canu-2.0.*.tar.xz -C tools
		"""

rule get_WTP:
	output:
		WTP_dir=directory(os.path.join(workflow.basedir, config['WTP_dir'])),
	message:
		"Downloading What the Phage"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_WTP/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p tools
		cd tools
		mkdir {output.WTP_dir}
		cd {output.WTP_dir}
		singularity pull  --name nanozoo-sourmash-3.4.1--16a8db7.img docker://nanozoo/sourmash:3.4.1--16a8db7
		nextflow run replikation/What_the_Phage -r v1.2.0 --setup -profile local,singularity --cachedir cache_dir
		"""

rule get_vcontact2:
	output:
		vcontact_dir=directory(config['vcontact_dir']),
	message:
		"Downloading vConTACT"
	conda:
		dirs_dict["ENVS_DIR"] + "/wtp.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/get_vcontact2/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.vcontact_dir}
		cd {output.vcontact_dir}
		wget https://bitbucket.org/MAVERICLab/vcontact2/src/master/vConTACT2.def
		singularity build --fakeroot vConTACT2.sif vConTACT2.def 
		"""

rule downloadDionSpacers:
	output:
		# blast=(os.path.join(workflow.basedir,"db/ncbi/NCBI_viral_proteins.faa")),
		dion_db=directory(os.path.join(workflow.basedir, config['dion_db'])),
	message:
		"Downloading Dion spacer database"
	threads: 1
	conda:
		dirs_dict["ENVS_DIR"]+ "/env1_spacepharer.yaml",
	benchmark:
		dirs_dict["BENCHMARKS"] + "/downloadDionSpacers/tot.tsv"
	shell:
		"""
		mkdir {output.dion_db}
		cd {output.dion_db}
		spacepharer downloaddb spacers_dion_et_al_2021 dionSetDB tmpFolder
		"""

rule download_bakta_db:
	output:
		db=directory(config["bakta_db"])
	message:
		"Downloading Bakta database"
	conda:
		dirs_dict["ENVS_DIR"] + "/bakta.yaml"
	benchmark:
		dirs_dict["BENCHMARKS"] + "/download_bakta_db/tot.tsv"
	threads: 1
	shell:
		"""
		mkdir -p {output.db}
		bakta_db download --output {output.db} 
		"""
