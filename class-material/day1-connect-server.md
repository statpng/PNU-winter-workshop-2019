---
layout: page
---

--- 

### II-1. Setting up the access to the compute server.

In this part, we will briefly learn how to access the KOBIC compute
server and to initially set up environments.

---

#### a. Make a connection to the server

- _IMPORTANT NOTE:_ Make sure that you have received the login information to access the KOBIC
  server.
- _IMPORTANT NOTE:_ For security reasons, we will denote the IP address of the
  server as `IP.AD.DR.ESS` and the login ID as `kobic00`. Please
  replace this information with your own login ID and IP address.
- The first thing you need to do is open a local terminal to access a
  shell prompt that you can type in.
  * If you have Windows laptop, open `MobaXTerm` and start a local
  terminal. Follow
  __[this link](https://mobaxterm.mobatek.net/download-home-edition.html)__ to
  download the software.
  * For Mac OS X users, just open the `Terminal` utility. If you do
  not know where it is located, open `Finder`, and click
  `Applications` on the left panel, and then open the `Utilities` folder,
  then you will be able to open `Terminal.app`. Keeping the
  `Terminal` application in the Dock would probably good for the next
  few days.
- To connect to the server, type `ssh kobic00@IP.AD.DR.ESS [Enter]` in
  your shell.
> <pre>
$ ssh kobic00@IP.AD.DR.ESS
kobic@@IP.AD.DR.ESS's password: [Need to type password here]
Last login: Sun Jul 15 05:59:10 2018 from IP.AD.DR.ESS
[kobic00@localhost ref]$  </pre>
  * Q. Were you able to successfully connect?
  * Q. What are the same and different from the example output shown
    above?

---

#### b. Make a symbolic link to the shared data directory

- Depending on users, the directory structures are slightly different.
  Notably, there are two replicated data directories, `/BiO/data`, and
  `/Bio2/data`. To avoid performance slowdown, it is important to
  use designated directories for each user.
- To ensure that each user access their designated data directory,
  we will create a symbolic link under your home directory to point
  the designated data directory.
  * If the last two digits of your username is between `01` and `25`,
  type the following commands. The second command is for verification
  purpose. Make sure that `BiO` is spelled correctly with lowercase
  `i`.
> <pre>
[kobic00@localhost ~]$ pwd
/BiO/home/kobic00
[kobic00@localhost ~]$ ln -s /BiO/data/ data
[kobic00@localhost ~]$ ls -l
total 0
lrwxrwxrwx 1 kobic00 kobic 11 Jul 15 08:08 data -> /BiO/data/ </pre>
  * If the last two digits of your username is between `26` and `50`,
  type the following commands. The second command is for verification
  purpose. Make sure that `BiO` is spelled correctly with lowercase `i`
> <pre>
[kobic00@localhost ~]$ pwd
/BiO2/home/kobic00
[kobic00@localhost ~]$ ln -s /BiO2/data/ data
[kobic00@localhost ~]$ ls -l
total 0
lrwxrwxrwx 1 kobic00 kobic 11 Jul 15 08:08 data -> /BiO2/data/ </pre>

---

Is everything OK? Otherwise, do not hesitate to ask questions to your
instructors, as this step is extremely important to be able to continue all
the next steps. 

If everything is okay, then let's go to next step [II-2 Processing raw sequence reads](../class-material/day1-fastq-practice.html)
, or go back to [Day 1 Overview](../day1).

---
