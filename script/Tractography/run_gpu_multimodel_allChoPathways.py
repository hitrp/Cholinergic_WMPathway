#!/usr/bin/env python3
"""Prepare HCP-style diffusion data and run the UKB NF Singularity pipeline."""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

import nibabel as nib
import numpy as np
import yaml


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_SIF = "data/unk_NF_multimodel_allChoPathways.sif"
CONTAINER_PYTHON = "tools/fsl/6.0.6.5/bin/python"
CONTAINER_PIPELINE = "/opt/pipeline/run_pipeline.py"
INTERNAL_TEMPLATES = "/opt/pipeline/resources/templates/crop_tone"
INTERNAL_CHOLINERGIC = "/opt/pipeline/resources/cholinergic"
REQUIRED_HCP_FILES = ("data.nii.gz", "bvals", "bvecs", "nodif_brain_mask.nii.gz")


class HcpRunnerError(RuntimeError):
    """A user-actionable HCP preparation or execution error."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare HCP-style {subject}/T1w/Diffusion inputs and run "
            "ukb_NF_unencrypted.sif."
        )
    )
    parser.add_argument(
        "--hcp-root",
        type=Path,
        required=True,
        help="Directory containing HCP subject directories.",
    )
    parser.add_argument(
        "--subjects",
        nargs="+",
        required=True,
        help="One or more HCP subject IDs.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        help="Pipeline output root (default: HCP root).",
    )
    parser.add_argument(
        "--template-masks-dir",
        type=Path,
        help="Override the SIF's embedded *_template.nii.gz masks.",
    )
    parser.add_argument(
        "--cholinergic-masks-dir",
        type=Path,
        help="Override the SIF's embedded cholinergic masks.",
    )
    parser.add_argument("--sif", type=Path, default=DEFAULT_SIF)
    parser.add_argument(
        "--singularity",
        default=None,
        help="Singularity/Apptainer executable (auto-detected by default).",
    )
    parser.add_argument(
        "--sudo",
        action="store_true",
        help="Run Singularity through sudo (only for hosts that require it).",
    )
    parser.add_argument("--gpu-device", type=int, default=0)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--stages", help="Comma-separated stages, for example 1,2.")
    selection.add_argument("--from-stage", type=int, choices=(1, 2, 3))
    parser.add_argument("--force", action="store_true", help="Rerun selected stages.")
    parser.add_argument(
        "--force-staging",
        action="store_true",
        help="Replace existing staged input links/files.",
    )
    parser.add_argument(
        "--prepare-only",
        action="store_true",
        help="Prepare inputs and configuration without running the SIF.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate and print the Singularity command without running it.",
    )
    parser.add_argument("--bval-max", type=float, default=1200.0)
    parser.add_argument("--bedpostx-fibres", type=int, default=3)
    parser.add_argument("--probtrack-samples", type=int, default=5000)
    return parser.parse_args()


def resolved_directory(path: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise HcpRunnerError(f"{label} does not exist or is not a directory: {resolved}")
    return resolved


def find_engine(explicit: str | None) -> str:
    if explicit:
        executable = shutil.which(explicit) if os.sep not in explicit else explicit
        if not executable or not Path(executable).is_file():
            raise HcpRunnerError(f"Singularity executable not found: {explicit}")
        return str(Path(executable).resolve())
    for name in ("apptainer", "singularity"):
        executable = shutil.which(name)
        if executable:
            return str(Path(executable).resolve())
    raise HcpRunnerError("Neither apptainer nor singularity is installed.")


def find_hcp_diffusion(hcp_root: Path, subject: str) -> Path:
    candidates = (
        hcp_root / subject / "T1w" / "Diffusion",
        hcp_root / subject / "Diffusion",
    )
    for candidate in candidates:
        if candidate.is_dir():
            missing = [name for name in REQUIRED_HCP_FILES if not (candidate / name).is_file()]
            if not missing:
                return candidate
    details = "; ".join(
        f"{path} (missing: "
        f"{', '.join(name for name in REQUIRED_HCP_FILES if not (path / name).is_file())})"
        for path in candidates
    )
    raise HcpRunnerError(f"HCP diffusion inputs are incomplete for {subject}: {details}")


def validate_hcp_images(source: Path) -> None:
    bvals = np.atleast_1d(np.loadtxt(source / "bvals", dtype=float))
    bvecs = np.atleast_2d(np.loadtxt(source / "bvecs", dtype=float))
    dwi = nib.load(str(source / "data.nii.gz"))
    mask = nib.load(str(source / "nodif_brain_mask.nii.gz"))

    if dwi.ndim != 4:
        raise HcpRunnerError(f"HCP DWI must be 4D: {source / 'data.nii.gz'} {dwi.shape}")
    if dwi.shape[3] != len(bvals):
        raise HcpRunnerError(
            f"HCP DWI/bvals mismatch in {source}: {dwi.shape[3]} vs {len(bvals)}"
        )
    if bvecs.shape != (3, len(bvals)):
        raise HcpRunnerError(
            f"HCP bvecs must have shape (3, {len(bvals)}) in {source}; "
            f"found {bvecs.shape}"
        )
    if mask.shape != dwi.shape[:3]:
        raise HcpRunnerError(
            f"HCP DWI/mask shape mismatch in {source}: {dwi.shape[:3]} vs {mask.shape}"
        )
    nodif = source / "nodif_brain.nii.gz"
    if nodif.is_file() and nib.load(str(nodif)).shape != dwi.shape[:3]:
        raise HcpRunnerError(f"HCP nodif brain shape does not match DWI: {nodif}")


def replace_symlink(destination: Path, target: str, force: bool) -> None:
    if destination.is_symlink():
        if os.readlink(destination) == target:
            return
        if not force:
            raise HcpRunnerError(
                f"Staged link points elsewhere: {destination}. Use --force-staging."
            )
        destination.unlink()
    elif destination.exists():
        if not force:
            raise HcpRunnerError(
                f"Staged input already exists: {destination}. Use --force-staging."
            )
        if destination.is_dir():
            raise HcpRunnerError(f"Refusing to replace directory: {destination}")
        destination.unlink()
    destination.symlink_to(target)


def generate_nodif_brain(source: Path, destination: Path, force: bool) -> None:
    if destination.is_file() and not force:
        return
    if destination.exists() or destination.is_symlink():
        if not force:
            raise HcpRunnerError(
                f"Generated nodif brain already exists: {destination}. "
                "Use --force-staging."
            )
        destination.unlink()

    bvals = np.atleast_1d(np.loadtxt(source / "bvals", dtype=float))
    b0_indices = np.flatnonzero(bvals <= 50)
    if not len(b0_indices):
        raise HcpRunnerError(f"No b=0 volumes (b <= 50) found in {source / 'bvals'}")

    dwi = nib.load(str(source / "data.nii.gz"))
    if dwi.ndim != 4 or dwi.shape[3] != len(bvals):
        raise HcpRunnerError(
            f"DWI/bvals mismatch for {source}: shape={dwi.shape}, bvals={len(bvals)}"
        )
    mask_image = nib.load(str(source / "nodif_brain_mask.nii.gz"))
    if mask_image.shape != dwi.shape[:3]:
        raise HcpRunnerError(
            f"DWI/mask shape mismatch for {source}: {dwi.shape[:3]} vs {mask_image.shape}"
        )

    mean_b0 = np.mean(np.asanyarray(dwi.dataobj)[..., b0_indices], axis=3)
    mask = np.asanyarray(mask_image.dataobj) > 0
    output = np.asarray(mean_b0 * mask, dtype=np.float32)
    header = dwi.header.copy()
    header.set_data_dtype(np.float32)
    nib.save(nib.Nifti1Image(output, dwi.affine, header), str(destination))


def stage_subject(
    hcp_root: Path,
    output_root: Path,
    subject: str,
    force: bool,
) -> tuple[Path, bool]:
    source = find_hcp_diffusion(hcp_root, subject)
    validate_hcp_images(source)
    destination = output_root / subject / "Diffusion" / "data"
    destination.mkdir(parents=True, exist_ok=True)
    relative_source = source.relative_to(hcp_root)
    container_source = Path("/hcp") / relative_source

    for name in REQUIRED_HCP_FILES:
        replace_symlink(destination / name, str(container_source / name), force)

    source_nodif = source / "nodif_brain.nii.gz"
    if source_nodif.is_file():
        replace_symlink(
            destination / "nodif_brain.nii.gz",
            str(container_source / "nodif_brain.nii.gz"),
            force,
        )
        generated = False
    else:
        generate_nodif_brain(source, destination / "nodif_brain.nii.gz", force)
        generated = True
    return source, generated


def validate_reference_files(templates: Path, cholinergic: Path) -> None:
    template_files = sorted(templates.glob("*_template.nii.gz"))
    if not template_files:
        raise HcpRunnerError(f"No *_template.nii.gz files found in {templates}")
    required_masks = ("AC.nii.gz", "Cingulum_ALL.nii.gz", "externalCapsule.nii.gz", "NbM.nii.gz")
    missing = [name for name in required_masks if not (cholinergic / name).is_file()]
    if missing:
        raise HcpRunnerError(
            f"Missing cholinergic masks in {cholinergic}: {', '.join(missing)}"
        )


def write_config(
    output_root: Path,
    subjects: list[str],
    gpu_device: int,
    args: argparse.Namespace,
    use_external_templates: bool,
    use_external_cholinergic: bool,
) -> Path:
    config_dir = output_root / ".ukb_nf"
    config_dir.mkdir(parents=True, exist_ok=True)
    config_file = config_dir / "hcp_pipeline_config.yaml"
    config = {
        "subjects": subjects,
        "paths": {
            "processed_root": "/output",
            "template_masks_dir": (
                "/refs/templates" if use_external_templates else INTERNAL_TEMPLATES
            ),
            "cholinergic_masks_dir": (
                "/refs/cholinergic"
                if use_external_cholinergic
                else INTERNAL_CHOLINERGIC
            ),
            "brainstem_mask": "/opt/qit/share/data/ants/masks/brainstem.nii.gz",
            "qit_atlas_dir": "/opt/qit/share/data",
            "atlas_fa_reference": "/opt/qit/share/data/crop/dti_FA.nii.gz",
            "logs_dir": "/output/logs",
        },
        "tools": {
            "qitdiff": "qitdiff",
            "qit": "qit",
            "flirt": "flirt",
            "fslmaths": "fslmaths",
            "bedpostx": "bedpostx_gpu",
            "probtrackx": "probtrackx2_gpu",
        },
        "execution": {"mock_external_tools": False, "gpu_device": gpu_device},
        "stage1": {
            "targets": [
                "diff.models.dti",
                "diff.models.fwdti",
                "diff.models.noddi",
                "diff.models.mcsmt",
                "atlas.models.dti",
                "atlas.models.fwdti",
                "atlas.models.noddi",
                "atlas.models.mcsmt",
            ],
            "checkpoints": {
                "diff.models.dti": "diff.models.dti/dti_FA.nii.gz",
                "diff.models.fwdti": "diff.models.fwdti/dti_FW.nii.gz",
                "diff.models.noddi": "diff.models.noddi/noddi_ficvf.nii.gz",
                "diff.models.mcsmt": "diff.models.mcsmt/mcsmt_diff.nii.gz",
                "atlas.models.dti": "atlas.models.dti/dti_FA.nii.gz",
                "atlas.models.fwdti": "atlas.models.fwdti/dti_FW.nii.gz",
                "atlas.models.noddi": "atlas.models.noddi/noddi_ficvf.nii.gz",
                "atlas.models.mcsmt": "atlas.models.mcsmt/mcsmt_diff.nii.gz",
            },
            "freesurfer_root": None,
            "extra_args": [],
        },
        "stage2": {
            "bval_max": args.bval_max,
            "bedpostx_fibres": args.bedpostx_fibres,
            "probtrack_samples": args.probtrack_samples,
        },
        "stage3": {
            "template_glob": "*_template.nii.gz",
            "registration_jvm_args": ["-Xmx5G"],
            "registration_seed": 42,
        },
    }
    temporary = config_file.with_suffix(f".tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    temporary.replace(config_file)
    return config_file


def build_command(
    engine: str,
    sif: Path,
    hcp_root: Path,
    output_root: Path,
    templates: Path | None,
    cholinergic: Path | None,
    args: argparse.Namespace,
) -> list[str]:
    command = []
    if args.sudo:
        sudo = shutil.which("sudo")
        if not sudo:
            raise HcpRunnerError("--sudo requested but sudo is not installed.")
        command.append(sudo)
    command.extend(
        [
            engine,
            "exec",
            "--cleanenv",
            "--nv",
            "--bind",
            f"{hcp_root}:/hcp:ro",
            "--bind",
            f"{output_root}:/output",
        ]
    )
    if templates is not None:
        command.extend(["--bind", f"{templates}:/refs/templates:ro"])
    if cholinergic is not None:
        command.extend(["--bind", f"{cholinergic}:/refs/cholinergic:ro"])
    command.extend(
        [
            str(sif),
            CONTAINER_PYTHON,
            CONTAINER_PIPELINE,
            "--config",
            "/output/.ukb_nf/hcp_pipeline_config.yaml",
        ]
    )
    if args.stages:
        command.extend(["--stages", args.stages])
    elif args.from_stage:
        command.extend(["--from-stage", str(args.from_stage)])
    if args.force:
        command.append("--force")
    return command


def main() -> int:
    args = parse_args()
    try:
        hcp_root = resolved_directory(args.hcp_root, "HCP root")
        output_root = (args.output_root or hcp_root).expanduser().resolve()
        output_root.mkdir(parents=True, exist_ok=True)
        templates = (
            resolved_directory(args.template_masks_dir, "Template masks directory")
            if args.template_masks_dir
            else None
        )
        cholinergic = (
            resolved_directory(args.cholinergic_masks_dir, "Cholinergic masks directory")
            if args.cholinergic_masks_dir
            else None
        )
        if templates is not None:
            template_files = sorted(templates.glob("*_template.nii.gz"))
            if not template_files:
                raise HcpRunnerError(f"No *_template.nii.gz files found in {templates}")
        if cholinergic is not None:
            required_masks = (
                "AC.nii.gz",
                "Cingulum_ALL.nii.gz",
                "externalCapsule.nii.gz",
                "NbM.nii.gz",
            )
            missing = [name for name in required_masks if not (cholinergic / name).is_file()]
            if missing:
                raise HcpRunnerError(
                    f"Missing cholinergic masks in {cholinergic}: {', '.join(missing)}"
                )

        sif = args.sif.expanduser().resolve()
        if not sif.is_file() or sif.stat().st_size == 0:
            raise HcpRunnerError(f"SIF does not exist or is empty: {sif}")
        if args.gpu_device < 0:
            raise HcpRunnerError("--gpu-device must be zero or greater.")
        if args.bedpostx_fibres < 1 or args.probtrack_samples < 1:
            raise HcpRunnerError("BEDPOSTX fibres and ProbtrackX samples must be positive.")

        subjects = list(dict.fromkeys(str(subject) for subject in args.subjects))
        print(f"HCP root: {hcp_root}")
        print(f"Output root: {output_root}")
        for subject in subjects:
            source, generated = stage_subject(
                hcp_root, output_root, subject, args.force_staging
            )
            suffix = "generated nodif_brain" if generated else "linked nodif_brain"
            print(f"Prepared {subject}: {source} ({suffix})")

        config_file = write_config(
            output_root,
            subjects,
            args.gpu_device,
            args,
            use_external_templates=templates is not None,
            use_external_cholinergic=cholinergic is not None,
        )
        print(f"Configuration: {config_file}")
        if args.prepare_only:
            print("Preparation complete; SIF execution skipped.")
            return 0

        engine = find_engine(args.singularity)
        command = build_command(
            engine, sif, hcp_root, output_root, templates, cholinergic, args
        )
        print(f"Command: {shlex.join(command)}")
        if args.dry_run:
            print("Dry run complete; SIF execution skipped.")
            return 0

        started = datetime.now()
        result = subprocess.run(command, check=False)
        elapsed = datetime.now() - started
        print(f"Singularity exit code: {result.returncode}; elapsed: {elapsed}")
        return result.returncode
    except (HcpRunnerError, OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
