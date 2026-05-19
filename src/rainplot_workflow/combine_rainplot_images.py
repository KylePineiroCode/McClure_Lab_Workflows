#!/usr/bin/env python3

import argparse
import os
from PIL import Image


def parse_args():
    """
    Parse command-line arguments for combining rain plot images with genomic
    annotation images.

    This script expects a manifest file that lists one or more rain plot PNGs.
    Each rain plot is stacked vertically with a genomic annotation PNG and saved
    into the output directory.

    The annotation image can either be:
    1. A shared annotation image passed through --annotation_png, or
    2. A matching annotation image located in the same directory as the rain plot.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments containing the manifest path, annotation
        PNG path, output directory, and whether input images should be deleted.
    """
    parser = argparse.ArgumentParser(
        description="Combine rain plot PNGs with a genomic feature PNG."
    )
    parser.add_argument(
        "--manifest",
        required=True,
        help="Path to a text file containing rain plot PNG paths, one per line."
    )
    parser.add_argument(
        "--annotation_png",
        required=True,
        help="Path to the genomic feature PNG."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        help="Directory to write combined PNGs to."
    )
    parser.add_argument(
        "--delete_inputs",
        action="store_true",
        help="Delete the source rain plot and annotation PNGs after combining."
    )
    return parser.parse_args()


def read_manifest(manifest_path):
    """
    Read rain plot PNG paths from a manifest file.

    The manifest file should contain one rain plot PNG path per line. Empty lines
    are ignored. After reading the paths, this function keeps only paths that
    actually exist on disk.

    This extra existence check is useful because a workflow may write a manifest
    before all expected plots are successfully created. Filtering missing files
    prevents the script from crashing later when trying to open a PNG that does
    not exist.

    Parameters
    ----------
    manifest_path : str
        Path to the text file containing rain plot PNG paths.

    Returns
    -------
    list[str]
        A list of existing rain plot PNG paths.
    """
    with open(manifest_path, "r", encoding="utf-8") as handle:
        # Strip whitespace/newlines and skip blank lines so the manifest can be
        # formatted cleanly without causing invalid empty file paths.
        paths = [line.strip() for line in handle if line.strip()]

    # Only return files that exist. This protects downstream image opening from
    # failing because of stale or missing paths in the manifest.
    return [path for path in paths if os.path.exists(path)]


def stack_images(top_path, bottom_path, output_path):
    """
    Stack two PNG images vertically and save the combined image.

    The first image is placed on top, and the second image is placed underneath.
    In this workflow, the top image is usually the rain plot and the bottom image
    is usually the genomic feature or annotation plot.

    The combined canvas uses the wider of the two images. Each image is centered
    horizontally on a white background so that images with different widths still
    align cleanly in the final output.

    Images are converted to RGB before combining to avoid issues with different
    PNG modes, such as transparency or palette-based images.

    Parameters
    ----------
    top_path : str
        Path to the image that should appear on the top of the combined figure.

    bottom_path : str
        Path to the image that should appear below the top image.

    output_path : str
        Path where the combined PNG should be written.
    """
    with Image.open(top_path) as top_image, Image.open(bottom_path) as bottom_image:
        # Convert both images to RGB so they share the same color mode.
        # This avoids errors when pasting/saving images that may have alpha
        # channels, palettes, or other PNG-specific modes.
        top = top_image.convert("RGB")
        bottom = bottom_image.convert("RGB")

        # The combined image must be wide enough for the wider of the two plots.
        # The height is the sum of both plots because we are stacking vertically.
        combined_width = max(top.width, bottom.width)
        combined_height = top.height + bottom.height

        # Use a white canvas because scientific plots usually have white
        # backgrounds, and this keeps any extra padding visually clean.
        canvas = Image.new("RGB", (combined_width, combined_height), "white")

        # Center each image horizontally. This matters when the rain plot and
        # annotation plot are not exactly the same width.
        canvas.paste(top, ((combined_width - top.width) // 2, 0))
        canvas.paste(bottom, ((combined_width - bottom.width) // 2, top.height))

        canvas.save(output_path)


def resolve_annotation_path(rainplot_path, shared_annotation_path):
    """
    Find the annotation image that should be paired with a rain plot.

    The script first checks whether there is a rain-plot-specific annotation PNG
    in the same directory as the rain plot. It does this by replacing the first
    occurrence of "rainplot_" in the rain plot filename with "annotation_".

    For example:
        rainplot_chrI_10000_20000.png

    would look for:
        annotation_chrI_10000_20000.png

    If that matching annotation image exists, it is used. Otherwise, the script
    falls back to the shared annotation image provided by --annotation_png.

    This approach supports both workflow styles:
    1. One annotation image shared across all rain plots.
    2. Separate annotation images for each rain plot/window.

    Parameters
    ----------
    rainplot_path : str
        Path to the current rain plot PNG.

    shared_annotation_path : str
        Path to the default annotation PNG provided by the user.

    Returns
    -------
    str
        Path to the annotation PNG that should be stacked with this rain plot.
    """
    annotation_name = os.path.basename(rainplot_path).replace(
        "rainplot_", "annotation_", 1
    )

    # Look for an annotation image in the same directory as the rain plot.
    # This allows each rain plot to have its own matching annotation track.
    candidate = os.path.join(os.path.dirname(rainplot_path), annotation_name)

    if os.path.exists(candidate):
        return candidate

    # If no rain-plot-specific annotation exists, use the shared annotation PNG.
    return shared_annotation_path


def main():
    """
    Run the full image-combining workflow.

    This function:
    1. Parses command-line arguments.
    2. Reads rain plot PNG paths from the manifest.
    3. Verifies that required input files exist.
    4. Creates the output directory if needed.
    5. Combines each rain plot with the correct annotation image.
    6. Optionally deletes the input PNGs after successful combining.

    The delete option is useful when this script is part of a larger workflow
    and the intermediate rain plot/annotation images are no longer needed after
    the final combined images are created.
    """
    args = parse_args()
    rainplot_paths = read_manifest(args.manifest)

    if not rainplot_paths:
        raise FileNotFoundError(
            f"No rain plot PNGs were found from manifest: {args.manifest}"
        )

    if not os.path.exists(args.annotation_png):
        raise FileNotFoundError(
            f"Genomic feature PNG not found: {args.annotation_png}"
        )

    # Create the output directory if it does not already exist.
    # exist_ok=True prevents an error if the directory is already there.
    os.makedirs(args.output_dir, exist_ok=True)

    for rainplot_path in rainplot_paths:
        rainplot_name = os.path.basename(rainplot_path)

        # Keep the original filename structure, but change the prefix so the
        # final image is clearly labeled as a combined plot.
        combined_name = rainplot_name.replace("rainplot_", "combined_", 1)
        output_path = os.path.join(args.output_dir, combined_name)

        # Use a specific annotation image if one matches this rain plot;
        # otherwise fall back to the shared annotation image.
        annotation_path = resolve_annotation_path(rainplot_path, args.annotation_png)

        stack_images(rainplot_path, annotation_path, output_path)
        print(f"[INFO] Combined plot written to: {output_path}", flush=True)

        # Delete the rain plot only after the combined image has been written.
        # This protects against losing input files if image creation fails.
        if args.delete_inputs and os.path.exists(rainplot_path):
            os.remove(rainplot_path)

        # Delete rain-plot-specific annotation images, but do not delete the
        # shared annotation here because it may still be needed for other plots.
        if (
            args.delete_inputs
            and annotation_path != args.annotation_png
            and os.path.exists(annotation_path)
        ):
            os.remove(annotation_path)

    # After all rain plots have been processed, it is safe to delete the shared
    # annotation image if the user requested cleanup.
    if args.delete_inputs and os.path.exists(args.annotation_png):
        os.remove(args.annotation_png)


if __name__ == "__main__":
    main()