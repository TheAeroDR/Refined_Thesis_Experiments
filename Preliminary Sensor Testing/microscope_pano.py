import cv2
import os
import numpy as np

def capture_on_keypress_and_create_pano(camera_index=0, output_folder="images/"):
    """
    Capture images on keypress and stitch them into a panorama.
    
    :param camera_index: Index of the camera to use.
    :param output_folder: Folder to save the captured images.
    """
    # Open the camera
    cap = cv2.VideoCapture(camera_index)

    if not cap.isOpened():
        print("Error: Could not open camera.")
        return
    
    # Create output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    images = []  # List to store captured images
    image_counter = 0

    while True:
        ret, frame = cap.read()
        if not ret:
            print("Error: Could not read frame.")
            break

        # Display the live feed
        cv2.imshow("Live Feed", frame)

        # Wait for a key press
        key = cv2.waitKey(1) & 0xFF

        # Capture the image when the 'c' key is pressed
        if key == ord('c'):
            image_counter += 1
            image_path = os.path.join(output_folder, f"captured_image_{image_counter}.jpg")
            cv2.imwrite(image_path, frame)
            images.append(frame)  # Add the captured frame to the images list
            print(f"Image captured and saved as {image_path}")

        # Exit the loop when 'q' key is pressed
        elif key == ord('q'):
            print("Exiting...")
            break

    # Release the camera and close all windows
    cap.release()
    cv2.destroyAllWindows()

def stitch_images_incrementally(images, visualize=False, mode='homography'):
    """
    Incrementally stitches a list of images into a panorama.
    
    :param images: List of images to stitch together.
    :param visualize: If True, show intermediate visualization windows for each step.
    :param mode: 'homography' (default) or 'translation' to force translation-only alignment.
    :return: Final stitched panorama image.
    """
    if len(images) < 2:
        print("At least 2 images are required to create a panorama.")
        return None
    
    print("Stitching images incrementally into a panorama...")
    pano = images[0]  # Start with the first image

    # Incrementally stitch images
    for i in range(1, len(images)):
        print(f"Stitching image {i+1} into panorama (step {i})...")
        pano = stitch_images(pano, images[i], visualize=visualize, step=i, mode=mode)
        if pano is None:
            print(f"Stitching failed between image {i} and {i+1}.")
            break

    return pano
def stitch_images(pano, image, visualize=False, step=0, mode='homography'):
    """
    Stitch a new image into the current panorama using feature matching + homography.
    Returns the updated panorama or None on failure.
    If visualize=True, shows intermediate images and waits for a keypress between steps.
    """
    # Ensure 3-channel BGR images
    if pano.ndim == 2:
        pano = cv2.cvtColor(pano, cv2.COLOR_GRAY2BGR)
    if image.ndim == 2:
        image = cv2.cvtColor(image, cv2.COLOR_GRAY2BGR)

    # Initialize SIFT detector
    sift = cv2.SIFT_create()

    # Detect keypoints and descriptors
    kp1, des1 = sift.detectAndCompute(pano, None)   # pano keypoints
    kp2, des2 = sift.detectAndCompute(image, None)  # new image keypoints

    if des1 is None or des2 is None:
        print("Not enough descriptors for matching.")
        return None

    # Use BFMatcher with kNN and Lowe's ratio test
    bf = cv2.BFMatcher()
    knn_matches = bf.knnMatch(des2, des1, k=2)  # match image -> pano
    good = []
    ratio = 0.75
    for m_n in knn_matches:
        if len(m_n) != 2:
            continue
        m, n = m_n
        if m.distance < ratio * n.distance:
            good.append(m)

    if len(good) < 4:
        print(f"Not enough good matches: {len(good)}")
        return None

    # Visualize keypoints and matches
    if visualize:
        kp_img = cv2.drawKeypoints(image, kp2, None, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
        kp_pano = cv2.drawKeypoints(pano, kp1, None, flags=cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
        cv2.imshow(f"Step {step}: Image keypoints", cv2.resize(kp_img, (min(kp_img.shape[1], 1200), min(kp_img.shape[0], 800))))
        cv2.imshow(f"Step {step}: Pano keypoints", cv2.resize(kp_pano, (min(kp_pano.shape[1], 1200), min(kp_pano.shape[0], 800))))
        matches_img = cv2.drawMatches(image, kp2, pano, kp1, good, None,
                                      flags=cv2.DrawMatchesFlags_NOT_DRAW_SINGLE_POINTS)
        cv2.imshow(f"Step {step}: Matches ({len(good)})", cv2.resize(matches_img, (min(matches_img.shape[1], 1400), min(matches_img.shape[0], 900))))
        k = cv2.waitKey(0) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
            return None
        cv2.destroyWindow(f"Step {step}: Image keypoints")
        cv2.destroyWindow(f"Step {step}: Pano keypoints")
        cv2.destroyWindow(f"Step {step}: Matches ({len(good)})")

    # Points from image (src) to pano (dst)
    src_pts = np.float32([kp2[m.queryIdx].pt for m in good]).reshape(-1, 1, 2)
    dst_pts = np.float32([kp1[m.trainIdx].pt for m in good]).reshape(-1, 1, 2)

    # choose mode
    if mode == 'translation':
        # Estimate pure translation by median offset of matched keypoints (dst - src).
        diffs = (dst_pts - src_pts).reshape(-1, 2)
        tx = float(np.median(diffs[:, 0]))
        ty = float(np.median(diffs[:, 1]))
        H = np.array([[1.0, 0.0, tx],
                      [0.0, 1.0, ty],
                      [0.0, 0.0, 1.0]], dtype=np.float64)
        # Plausibility check
        max_dim = max(pano.shape[1], pano.shape[0], image.shape[1], image.shape[0])
        if abs(tx) > max_dim * 2 or abs(ty) > max_dim * 2:
            print(f"Estimated translation implausible (tx={tx:.1f}, ty={ty:.1f}), aborting.")
            return None
    else:
        # Compute homography that maps image -> pano
        H, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
        if H is None:
            print("Homography computation failed.")
            return None

    h1, w1 = pano.shape[:2]
    h2, w2 = image.shape[:2]

    # Compute transformed corners to determine output canvas size
    corners_img = np.float32([[0, 0], [w2, 0], [w2, h2], [0, h2]]).reshape(-1, 1, 2)
    warped_corners = cv2.perspectiveTransform(corners_img, H)
    pano_corners = np.float32([[0, 0], [w1, 0], [w1, h1], [0, h1]]).reshape(-1, 1, 2)

    all_corners = np.concatenate((warped_corners, pano_corners), axis=0)
    [xmin, ymin] = np.int32(all_corners.min(axis=0).ravel() - 0.5)
    [xmax, ymax] = np.int32(all_corners.max(axis=0).ravel() + 0.5)

    translate = [-xmin, -ymin]
    H_trans = np.array([[1, 0, translate[0]],
                        [0, 1, translate[1]],
                        [0, 0, 1]]) @ H

    out_width = xmax - xmin
    out_height = ymax - ymin

    # Warp the new image into the panorama coordinate frame
    warped_image = cv2.warpPerspective(image, H_trans, (out_width, out_height))

    # Visualize warped image and transformed corners
    if visualize:
        cv2.imshow(f"Step {step}: Warped image", cv2.resize(warped_image, (min(warped_image.shape[1], 1400), min(warped_image.shape[0], 900))))
        k = cv2.waitKey(0) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
            return None
        cv2.destroyWindow(f"Step {step}: Warped image")

    # Create canvas and place pano into it
    canvas = np.zeros((out_height, out_width, 3), dtype=np.uint8)
    x_off, y_off = translate
    canvas[y_off:y_off + h1, x_off:x_off + w1] = pano

    # Visualize canvas before blending
    if visualize:
        cv2.imshow(f"Step {step}: Canvas before blend", cv2.resize(canvas, (min(canvas.shape[1], 1400), min(canvas.shape[0], 900))))
        k = cv2.waitKey(0) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
            return None
        cv2.destroyWindow(f"Step {step}: Canvas before blend")

        # Create masks and blend where both images contribute
    mask_canvas = (canvas > 0).any(axis=2)
    mask_warp = (warped_image > 0).any(axis=2)

    overlap = mask_canvas & mask_warp
    only_canvas = mask_canvas & ~mask_warp
    only_warp = mask_warp & ~mask_canvas

    # Visualize overlap mask
    if visualize:
        mask_vis = np.zeros_like(canvas)
        mask_vis[only_canvas] = (0, 255, 0)   # green = only pano
        mask_vis[only_warp] = (255, 0, 0)     # blue = only warped image
        mask_vis[overlap] = (0, 0, 255)       # red = overlap
        cv2.imshow(f"Step {step}: Contribution mask", cv2.resize(mask_vis, (min(mask_vis.shape[1], 1400), min(mask_vis.shape[0], 900))))
        k = cv2.waitKey(0) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
            return None
        cv2.destroyWindow(f"Step {step}: Contribution mask")

    # Convert to float for blending
    canvas_f = canvas.astype(np.float32)
    warp_f = warped_image.astype(np.float32)

    # --- Improved blending: exposure compensation + distance-transform feathering ---
    # Simple exposure compensation: match mean color in the overlap
    if np.any(overlap):
        eps = 1e-8
        mean_canvas = np.array([canvas_f[..., c][overlap].mean() if np.any(overlap) else 1.0 for c in range(3)])
        mean_warp = np.array([warp_f[..., c][overlap].mean() if np.any(overlap) else 1.0 for c in range(3)])
        gain = mean_canvas / (mean_warp + eps)
        # Limit gain to avoid extreme scaling
        gain = np.clip(gain, 0.7, 1.3)
        warp_f = warp_f * gain.reshape(1, 1, 3)

    # Create smooth alpha using distance transform so seam fades gradually
    mask_canvas_u8 = (mask_canvas.astype(np.uint8) * 255)
    mask_warp_u8 = (mask_warp.astype(np.uint8) * 255)

    # Distance transform (zero outside the mask)
    dist_canvas = cv2.distanceTransform(mask_canvas_u8, cv2.DIST_L2, 5).astype(np.float32)
    dist_warp = cv2.distanceTransform(mask_warp_u8, cv2.DIST_L2, 5).astype(np.float32)

    alpha_warp = np.zeros((out_height, out_width), dtype=np.float32)
    sum_dist = dist_warp + dist_canvas
    # Where both present, weight by relative distance to boundary
    overlap_idx = overlap
    alpha_warp[overlap_idx] = dist_warp[overlap_idx] / (sum_dist[overlap_idx] + 1e-8)
    # Only-warp and only-canvas regions
    alpha_warp[only_warp] = 1.0
    alpha_warp[only_canvas] = 0.0

    # Expand alpha to 3 channels and blend
    alpha_3 = np.repeat(alpha_warp[:, :, np.newaxis], 3, axis=2)
    result_f = canvas_f * (1.0 - alpha_3) + warp_f * alpha_3

    result = np.clip(result_f, 0, 255).astype(np.uint8)

    # Visualize final stitched result for this step
    if visualize:
        cv2.imshow(f"Step {step}: Result", cv2.resize(result, (min(result.shape[1], 1400), min(result.shape[0], 900))))
        k = cv2.waitKey(0) & 0xFF
        if k == ord('q'):
            cv2.destroyAllWindows()
            return None
        cv2.destroyWindow(f"Step {step}: Result")

    return result

def main():
    # Call the function to start live feed and capture images on keypress
    capture_on_keypress_and_create_pano(camera_index=12, output_folder="images/")
    images = []

    # os hunt for 'captured_image_X.jpg' in images/ folder
    # need to specially load them in ascending numeric order
    for root, dirs, files in os.walk("images/"):
        files.sort(key=lambda f: int(os.path.splitext(f)[0].split('_')[-1]))
        for file in files:
            if file.startswith("captured_image_") and file.endswith(".jpg"):
                images.append(cv2.imread(os.path.join(root, file)))

    if images:
        stitch_mode = 'translation'  # 'homography' or 'translation'
        panorama = stitch_images_incrementally(images, visualize=False, mode=stitch_mode)
        if panorama is not None:
            cv2.imshow("Final Panorama", cv2.resize(panorama, (min(panorama.shape[1], 1400), min(panorama.shape[0], 900))))
            #wait for 5s
            cv2.waitKey(5000)
            cv2.destroyAllWindows()
            cv2.imwrite("microscope_panorama.jpg", panorama)
            print("Panorama created and saved as microscope_panorama.jpg")
        else:
            print("Failed to create panorama.")
    else:
        print("No images to stitch.")

if __name__ == "__main__":
    main()
