# Homework 3 Grading Script
We will use this script to grade your program. **Make sure your program can be executed by this script.**

## Preparing
* Step1:  
    Enter the `HW3_grading` directory and create a new directory named with your student ID in the `student` directory.
    ```sh
    $ cd HW3_grading/
    $ mkdir student/${your_student_id}
    ```
    For example,
    ```sh
    $ cd HW3_grading/
    $ mkdir student/113000000
    ```

* Step2:  
    Put your compressed file in the directory which you just created.  
    The correct path should be:
    ```
    HW3_grading/student/${your_student_id}/CS6135_HW3_${your_student_id}.tar.gz
    ```
    For example,
    ```
    HW3_grading/student/113000000/CS6135_HW3_113000000.tar.gz
    ```

### Notice:  
**Please make sure not to put your original directory here**, as it will remove all directories before unzipping the compressed file.

## Grading
* Step1:  
    Navigate to the `HW3_grading` directory and run `HW3_grading.sh`.
    ```sh
    $ cd HW3_grading/
    $ bash HW3_grading.sh
    ```

* Step3:  
    Check your output. Ensure the status of each checking item is **yes**.
    * If the status of a testcase is **success**, it means your program finished in time, and the output result is legal.
        ```
        host name: ic51
        compiler version: g++ (GCC) 9.3.0

        grading on 113000000:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes

          testcase |      ratio | wirelength |    runtime | status
        -----------|------------|------------|------------|--------
           public1 |       0.15 |     199767 |       1.66 | success
           public2 |       0.15 |     388850 |       6.82 | success
           public3 |       0.15 |     517428 |      14.24 | success
           public1 |        0.1 |     225666 |       3.76 | success
           public2 |        0.1 |     407635 |      10.27 | success
           public3 |        0.1 |     560473 |      20.91 | success
        ```

    * If the status of a testcase is not **success**, it means your program failed in this testcase.
        ```
        host name: ic51
        compiler version: g++ (GCC) 9.3.0
        
        grading on 113000000:
         checking item          | status
        ------------------------|--------
         correct tar.gz         | yes
         correct file structure | yes
         have README            | yes
         have Makefile          | yes
         correct make clean     | yes
         correct make           | yes
        
          testcase |      ratio | wirelength |    runtime | status
        -----------|------------|------------|------------|--------
           public1 |       0.15 |        N/A |        TLE | Time out while testing public1.
           public2 |       0.15 |        N/A |        TLE | Time out while testing public2.
           public3 |       0.15 |        N/A |        TLE | Time out while testing public3.
           public1 |        0.1 |        N/A |        TLE | Time out while testing public1.
           public2 |        0.1 |        N/A |        TLE | Time out while testing public2.
           public3 |        0.1 |        N/A |        TLE | Time out while testing public3.
        ```
