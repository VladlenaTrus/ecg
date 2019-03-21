#include <stdio.h>
#include "stdlib.h"
#include <string.h>


unsigned short get16(char *p)
{
    unsigned short rez;
    rez  = *p++ & 0xff;
    rez |= *p << 8;
    return (rez);
}


unsigned long get32(char *p)
{
    unsigned long rez;
    rez  =  *p++ & 0xff;
    rez |= (*p++ & 0xff) << 8;
    rez |= (*p++ & 0xff) << 16;
    rez |=  *p << 24;
    return (rez);
}

int flag_start = 0;


char *name_of_prog;       /* the name of this program (for use in error messages) */
int num_of_leads = -1;    /* number of ECG leads */
int nqrs = -1;        /* number of beats detected */
unsigned long *nsamp;    /* nsamp[i]: number of samples in lead i */
long rblenms;        /* length of refbeat in milliseconds */
int fcM;        /* index of fiducial within refbeat */
char date_of_rec[20];    /* date of recording */
char time_of_record[20];    /* time of recording */
int blfilt = -1;    /* baseline filter -3 dB point */
char *inst = "<missing>";/* institution description */
int lpfilt = -1;    /* low-pass filter -3 dB point */
int filter_bit = -1;    /* see SCP spec */
char nm_patient[100];    /* patient's name (LAST, FIRST) */
char referring_dr[100];    /* referring dr */
char comments[200];    /* comments */
char id_patient[30];    /* patient ID (medical record number) */
int age = -1;        /* patient's age */
int sex = -1;        /* 1: male, 2: female, 0: unknown, 9: unspecified */
int RR_interval = -1;         /* RR interval */
int tdevid = -1;    /* ID number of acquiring device */
int tinst = -1;         /* institution number */
int tdept = -1;         /* department number */
struct subzone {    /* reference beat subtraction zone */
    short type;        /* beat type (0: dominant) */
    long t0;        /* first sample in zone */
    long t1;        /* fiducial */
    long t2;        /* last sample in zone */
} *subz;

short **ecg;        /* ecg[i]: reconstructed ECG for ith lead */
short **refbeat;    /* reference beat (template) */




main(int argc, char **argv)
{
    unsigned short crc, id_of_sec;
    unsigned long btsread, length, sec_len;
    unsigned char desc[20], header[6], *data, *p;
    int i;
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 's':
                    flag_start++;
                    break;
                default:
                    fprintf(stderr, "%s: unrecognized option '%s'\n", name_of_prog, argv[i]);
                    exit(1);
            }
        } else {
            fprintf(stderr,"%s: unrecognized arg '%s'\n", name_of_prog, argv[i]);
            exit(1);
        }
    }

    btsread = fread(header, 1, 6, stdin);
    if (btsread != 6) {
        fprintf(stderr, "%s: input is too short (%lu bytes)\n",
                name_of_prog, btsread);
        exit(1);
    }
    crc = get16(header);
    length = get32(header+2);
    if (flag_start) {
        printf("CRC = %d, expected length = %ld bytes\n", crc, length);
    }
    if ((data = (unsigned char *)malloc(length)) == NULL) {
        fprintf(stderr, "%s: not enough memory (%ld bytes needed)\n",
                name_of_prog, length);
        exit(2);
    }
    memcpy(data, header, 6);
    i = 1;
    while (!ferror(stdin) && btsread < length && i > 0) {
        btsread += i = fread(data + btsread, 1, length - btsread, stdin);
    }
    if (btsread < length) {
        fprintf(stderr, "%s: input too short (%lu byte%s missing)\n",
                name_of_prog, length-btsread, (length-btsread == 1) ? "" : "s");
        exit(3);
    }
    if (flag_start) {
        printf("Record length: %ld bytes\n", length);
    }

    p = data + 6;

    while (p < data + length) {
        crc = get16(p);
        id_of_sec = get16(p+2);
        sec_len = get32(p+4);
        if (flag_start) {
            printf("\nSection %d length: %ld bytes\n", id_of_sec, sec_len);
        }
        if (sec_len < 8) {
            if (flag_start) printf(
                        " Warning: section length too short (must be at least 8 bytes)\n"
                        "  Remaining data (%lu bytes) following short section will not be read\n",
                        (unsigned long) (data+length-(p+8)));
            sec_len = 8;
        }
        if (sec_len > data + length - p) {
            if (flag_start) {
                printf(

                        " Warning: section length exceeds amount of remaining data\n"
                        "  This section will not be parsed\n");
            }
            break;
        }
        sprintf(desc, " Section %d CRC", id_of_sec);
        switch (id_of_sec) {
            case 0: if (flag_start) printf("Contents:  Pointers to data areas in the record\n");
                section0(p, sec_len);
                break;
            case 1: if (flag_start) printf("Contents:  Header information - Patient data/ECG acquisition"" data\n");
                section1(p, sec_len);
                break;
            case 4: if (flag_start) printf("Contents:  QRS locations\n");
                section4(p, sec_len);
                break;
            case 7: if (flag_start) printf("Contents:  Global measurements\n");
                section7(p, sec_len);
                break;
            default:
                break;
        }
        p += sec_len;
    }


    if (refbeat) {
        for (i = 0; i < num_of_leads && refbeat[i] != NULL; i++) {
            free(refbeat[i]);
        }
        free(refbeat);
    }
    if (ecg) {
        for (i = 0; i <= num_of_leads && ecg[i] != NULL; i++) {
            free(ecg[i]);
        }
        free(ecg);
    }
    if (nsamp) {
        free(nsamp);
    }
    if (subz) {
        free(subz);
    }
    if (data) {
        free(data);
    }

    exit(0);
}


int section0(unsigned char *p, long len)
{
    p += 16; len -= 16;

    while (len > 0) {
        unsigned short sec_id;
        unsigned long sec_len, sec_index;

        sec_id = get16(p);
        sec_len = get32(p+2);
        sec_index = get32(p+6);
        if (flag_start) {
            printf("  Section %2d:%6ld bytes beginning at byte %6ld\n",
                   sec_id, sec_len, sec_index);
        }
        len -= 10;
        p += 10;
    }
    return (1);
}

char *month[] = { "January", "February", "March", "April", "May", "June",
                  "July", "August", "September", "October", "November", "December" };


int section1(unsigned char *p, long len) {
    FILE *pfile;
    unsigned short vlen;
    char format[16], *tpatid = NULL, *trefdr = NULL, *tcomments = NULL,
            *tdate = NULL, *ttime = NULL, *pyoa = NULL, *pyob = NULL;
    char *firstname = NULL, *lastname = NULL;
    enum { text, date, time, b1, b2, b2b1, bn, mixed } vtype;
    unsigned char *p_save = p;
    long len_save = len;

    p += 16; len -= 16;


    while (len > 2) {
        if (*p == 255) {
            break;
        }
        vlen = get16(p+1);
        switch (*p) {
            case  0:
                if (flag_start) {
                    printf("\n  Last name\n");
                }
                vtype = text; lastname = p+3;
                break;
            case  1:
                if (flag_start) {
                printf("\n  First name\n");
            }
                vtype = text; firstname = p+3;
                break;
            case  2:
                if (flag_start) {
                    printf("\n  Patient ID number\n");
                }
                vtype = text; tpatid = p+3;
                break;
            case  3:
                if (flag_start) {
                    printf("\n  Second last name\n");
                }
                vtype = text;
                break;
            case  4:
                if (flag_start) {
                    printf("\n  Age\n");
                }
                vtype = b2b1; age = get16(p+3);
                break;
            case  5:
                if (flag_start) {
                    printf("\n  Date of birth\n");
                }
                vtype = date;
                break;
            case  6:
                if (flag_start) {
                    printf("\n  Height\n");
                }
                vtype = b2b1; break;
            case  7:
                if (flag_start) {
                    printf("\n  Weight\n");
                }
                vtype = b2b1; break;
            case  8:
                if (flag_start) {
                    printf("\n  Sex\n");
                }
                vtype = b1; sex = *(p+3); break;
            case  9:
                if (flag_start) {
                    printf("\n  Race\n");
                }
                vtype = b1; break;
            case 10:
                if (flag_start) {
                    printf("\n  Drug\n");
                }
                vtype = mixed; break;
            case 11:
                if (flag_start) {
                    printf("\n  Systolic blood pressure\n");
                }
                vtype = b2; break;
            case 12:
                if (flag_start) {
                    printf("\n  Diastolic blood pressure\n");
                }
                vtype = b2; break;
            case 13:
                if (flag_start) {
                    printf("\n  Diagnosis or referral indication\n");
                }
                vtype = text; break;
            case 14:
                if (flag_start) {
                    printf("\n  ID of the acquiring device\n");
                }
                tdevid = get16(p+7);
                tinst = get16(p+3);
                tdept = get16(p+5);
                if (flag_start) {
                    char *q;
                    printf("   Institution number: %d\n", get16(p+3));
                    printf("   Department number: %d\n", get16(p+5));
                    printf("   Device ID: %d\n", get16(p+7));
                    printf("   Device type: %d (%s)\n", *(p+9),
                           (*(p+9) == 0 ? "cart" :
                            (*(p+9) == 1 ? "system or host" : "<error>")));
                    printf("   Manufacturer code: %d\n", *(p+10));
                    printf("   Text model description: [%d]", *(p+11));
                    printf("   SCP-ECG protocol revision number: %g\n",
                           (double)(*(p+17))/10.0);
                    printf("   SCP-ECG protocol compatibility level: %d\n",
                           *(p+18));
                    printf("   Language support code: %d\n", *(p+19));
                    printf("   Device capabilities: %d\n", *(p+20));
                    printf("   AC mains frequency environment: %d ", *(p+21));
                    switch (*(p+21)) {
                        case 0: printf("(unspecified)\n"); break;
                        case 1: printf("(50 Hz)\n"); break;
                        case 2: printf("(60 Hz)\n"); break;
                        default: printf("<error>"); break;
                    }
                    q = p+39;
                    printf("   Analyzing program revision number: [%s]\n", q);
                    q += strlen(q) + 1;
                    printf("   Serial number of acquisition device: [%s]\n", q);
                    q += strlen(q) + 1;
                    printf("   Acquisition device system software ID: [%s]\n", q);
                    q += strlen(q) + 1;
                    printf("   Acquisition device SCP implementation ID: [%s]\n",
                           q);
                    q += strlen(q) + 1;
                    printf("   Manufacturer of acquisition device: [%s]\n", q);
                }
                vtype = mixed; break;
            case 15:
                if (flag_start) {
                    printf("\n  ID of the analyzing device\n");
                }
                vtype = mixed; break;
            case 16:
                if (flag_start) {
                    printf("\n  Acquiring institution\n");
                }
                vtype = text; inst=p+3;
                break;
            case 17:
                if (flag_start) {
                    printf("\n  Analyzing institution\n");
                }
                vtype = text;
                break;
            case 18:
                if (flag_start) {
                    printf("\n  Acquiring department\n");
                }
                vtype = text;
            case 19:
                if (flag_start) {
                    printf("\n  Analyzing department\n");
                }
                vtype = text;
                break;
            case 20: if (flag_start) printf("\n  Referring physician\n");
                vtype = text;
                trefdr = p+3; break;
            case 21: if (flag_start) printf("\n  Latest confirming physician\n");
                vtype = text;
                break;
            case 22: if (flag_start) printf("\n  Technician\n");
                vtype = text;
                break;
            case 23: if (flag_start) printf("\n  Room\n");
                vtype = text;
                break;
            case 24: if (flag_start) printf("\n  Stat code (urgency)\n");
                vtype = b1; break;
            case 25: if (flag_start) printf("\n  Date of acquisition\n"); tdate = p+3;
                vtype = date; break;
            case 26: if (flag_start) printf("\n  Time of acquisition\n"); ttime = p+3;
                vtype = time; break;
            case 27: if (flag_start) printf("\n  Baseline filter\n");
                vtype = b2; blfilt = get16(p+3); break;
            case 28: if (flag_start) printf("\n  Low-pass filter\n");
                vtype = b2; lpfilt = get16(p+3); break;
            case 29: if (flag_start) printf("\n  Filter bit map\n");
                vtype = b1; filter_bit = *(p+3); break;
            case 30: if (flag_start) printf("\n  Free text field (comments)\n"); tcomments = p+3;
                vtype = text; break;
            case 31: if (flag_start) printf("\n  Sequence number\n");
                vtype = text; break;
            case 32: if (flag_start) printf("\n  Medical history codes\n");
                vtype = bn; break;
            case 33: if (flag_start) printf("\n  Electrode configuration code\n");
                vtype = b2; break;
            case 34: if (flag_start) printf("\n  Time zone\n");
                vtype = mixed; break;
            case 35: if (flag_start) printf("\n  Free text medical history\n");
                vtype = text; break;
            default: if (flag_start) printf("\n  <Undefined tag>\n");
                vtype = mixed; break;
        }
        if (flag_start) printf("   Tag%3d (%2d bytes): ", *p, vlen);
        p += 3; len -= vlen+3;
        if (vlen == 0) {
            if (flag_start) {
                printf("<not defined>\n");
            }
        }
        else switch (vtype) {
                case text:
                    if (vlen > 0) {
                        sprintf(format, "[%%%ds]", vlen-1);
                        if (flag_start) printf(format, p);
                        p += vlen;
                    }
                    else
                    if (flag_start) printf("<undefined>");
                    break;
                case date:
                    if (vlen != 4) {
                        if (flag_start) printf("  Error: incorrect date format (%d bytes)\n", vlen);
                    }
                    else {
                        if (flag_start) printf("%4d/%02d/%02d", get16(p), *(p+2), *(p+3));
                    }
                    p += vlen;
                    break;
                case time:
                    if (vlen != 3) {
                        if (flag_start) printf("  Error: incorrect time format (%d bytes)\n", vlen);
                    }
                    else
                    if (flag_start) printf("%02d:%02d:%02d", *p, *(p+1), *(p+2));
                    p += vlen;
                    break;
                case b1:
                    if (vlen != 1) {
                        if (flag_start) printf("  Error: incorrect format (%d bytes instead of 1)\n",
                                          vlen);
                    }
                    else
                    if (flag_start) printf("%d", *p);
                    p += vlen;
                    break;
                case b2:
                    if (vlen != 2) {
                        if (flag_start) printf("  Error: incorrect format (%d bytes instead of 2)\n",
                                          vlen);
                    }
                    else
                    if (flag_start) printf("%d", get16(p));
                    p += vlen;
                    break;
                case b2b1:
                    if (vlen != 3) {
                        if (flag_start) printf("  Error: incorrect format (%d bytes instead of 3)\n",
                                          vlen);
                    }
                    else
                    if (flag_start) printf("%d (code %d)", get16(p), *(p+2));
                    p += vlen;
                    break;
                case bn:
                    while (vlen-- > 0) {
                        if (flag_start) printf("%d, ", *p);
                        p++;
                    }
                    break;
                default:
                    if (flag_start) printf("oops!  undefined variable type\n");
                case mixed:
                    if (flag_start) {
                        while (vlen-- > 0) {
                            if (isprint(*p)) putchar(*p);
                            else if (flag_start) printf("\\%03o", *p);
                            p++;
                        }
                        printf("\n");
                    }
                    else
                        p += vlen;
            }
    }

    sprintf(date_of_rec, "%d %s %d",*(tdate+3), month[*(tdate+2)-1], get16(tdate));
    sprintf(time_of_record, "%02d:%02d:%02d", *ttime, *(ttime+1), *(ttime+2));
    if (lastname == NULL) lastname = "<last name missing>";
    if (firstname == NULL) firstname = "<first name missing>";
    if (strlen(lastname) + strlen(firstname) < sizeof(nm_patient) - 3)
        sprintf(nm_patient, "%s, %s", lastname, firstname);
    else {
        fprintf(stderr, "Patient name (%s, %s) is too long\n", lastname,
                firstname);
        strcpy(nm_patient, "Name not available -- too long");
    }
    if (tpatid == NULL) tpatid = "<patient ID missing>";
    strncpy(id_patient, tpatid, sizeof(id_patient));
    if (trefdr == NULL) trefdr = "<referring physician missing>";
    strncpy(referring_dr, trefdr, sizeof(referring_dr));
    if (tcomments == NULL) tcomments = "<comments missing>";
    strncpy(comments, tcomments, sizeof(comments));

    if (*p != 255)
        if (flag_start) printf(" Error: header terminator (tag 255) is missing\n");
    if (tpatid == NULL)
        if (flag_start) printf(" Error: patient ID number (tag 2) is missing\n");
    if (tdevid == -1 || tinst == -1 || tdept == -1)
        if (flag_start) printf(" Error: device ID number (tag 14) is missing\n");
    if (tdate == NULL)
        if (flag_start) printf(" Error: date of acquisition (tag 25) is missing\n");


    if (*p != 255 || tpatid == NULL || tdevid == -1 ||
        tdate == NULL || ttime == NULL) {
        return (0);
    }
    return (1);
}



int section4(unsigned char *p, long len)
{
    unsigned short i, stat = 1;

    p += 16; len -= 16; /* move to data area */
    if (len < 6) {
        if (flag_start) printf("  Error: section 4 contains no data\n");
        return (0);
    }
    rblenms = get16(p);
    fcM = get16(p+2);
    nqrs = get16(p+4);


    if (flag_start) {
        printf("  Length of reference beat type 0 data: %ld msec\n", rblenms);
        printf(
                "  Sample # of fiducial relative to start of ref beat type 0: %d\n",
                fcM);
        printf("  Number of QRS complexes in the entire record: %d\n", nqrs);
    }

    if (nqrs < 1)
        return (0);    /* no beats -- nothing else to do here */

    p += 6; len -= 6;
    if ((subz=(struct subzone *)malloc(sizeof(struct subzone)*nqrs)) == NULL) {
        fprintf(stderr, "  Error: too many (%d) beats\n", nqrs);
        return (0);
    }

    if (len >= nqrs*14) {
        if (flag_start) {
            printf("  Reference beat subtraction zones:\n");
        }
        for (i = 0; i < nqrs && len >= 14; i++) {
            subz[i].type = get16(p);
            subz[i].t0 = get32(p+2);
            subz[i].t1 = get32(p+6);
            subz[i].t2 = get32(p+10);
            if ((subz[i].type == 0) && (subz[i].t2 == 0)) {
                subz[i].type = 9;
            }

            if (flag_start) {
                printf(
                        "   %2d  Type: %2d  start: %7ld"
                        "  fiducial: %7ld  end: %7ld\n",
                        i+1, subz[i].type, subz[i].t0, subz[i].t1, subz[i].t2);
                if (subz[i].type != 0 && (subz[i].t0 != 0 ||
                                          subz[i].t2 != 0))
                    printf("  Error: start and end should be zero for beats "
                           "of non-zero types\n");
            }
            p += 14; len -= 14;
        }
    }
    else {
        if (flag_start) {
            printf("  Error: reference beat subtraction zones are missing\n");
        }
        stat = 0;
    }
    if (len >= nqrs * 8) {
        if (flag_start) {
            printf("  QRS locations:\n");
        }
        for (i = 0; i < nqrs; i++) {
            if (flag_start) {
                printf("   %2d  start: %7lu  end: %7lu\n",
                       i + 1, get32(p), get32(p + 4));
            }
            p += 8; len -= 8;
        }
    }
    else {
        if (flag_start) {
            printf("  Error: QRS locations are missing\n");
        }
        stat = 0;
    }
    return (stat);
}



int section7(unsigned char *p, long len)
{
    unsigned char nmb, nps, sec7qrs, *q, *r;
    short i, n, sec7qrs_offset;

    p += 16; len -= 16; /* move to data area */


    RR_interval = get16(p+2);


    nmb = *p;
    nps = *(p+1);
    sec7qrs_offset = 6 + 16 * nmb + 10 * nps;
    if (sec7qrs_offset > len) {
        sec7qrs_offset = len;
    }

    sec7qrs = get16(p + sec7qrs_offset);

    if (flag_start) {
        printf(" QRS complexes measured: %d\n", sec7qrs);
        printf(" Measurement blocks: %d\n", nmb);
        if (nmb == sec7qrs + 1)
            printf("  (First measurement block contains measurements for"
                   " reference beat type 0;\n"
                   "   Others contain measurements for each individual beat)\n");
        else
            printf("  (Measurements for each reference beat type)\n");
        printf("  Number of measurement blocks: %d\n", nmb);
        printf("  Number of pacemaker spikes: %d\n", nps);
        printf("  Mean RR interval: %d ms\n", RR_interval);
        printf("  Mean PP interval: %d ms\n", get16(p+4));
        for (i = 0, q = p+6; i < nmb; i++, q += 16) {
            printf("   Block %d\n", i);
            printf("    P onset:   %d ms\n", get16(q));
            printf("    P end:     %d ms\n", get16(q+2));
            printf("    QRS onset: %d ms\n", get16(q+4));
            printf("    QRS  end:  %d ms\n", get16(q+6));
            printf("    T end:     %d ms\n", get16(q+8));
            if ((n = get16(q+10)) == 999) printf("    P axis:    undefined\n");
            else printf("    P axis:    %d degrees\n", n);
            if ((n = get16(q+12)) == 999) printf("    QRS axis:  undefined\n");
            else printf("    QRS axis:  %d degrees\n", n);
            if ((n = get16(q+14)) == 999) printf("    T axis:    undefined\n");
            printf("    T axis:    %d degrees\n", n);
        }


        printf(" sec7qrs = <%d>\n", sec7qrs);
        if (sec7qrs > 0)printf(" QRS type information (%d beats)\n", get16(r));
        for (i = 0, q = r+2; i < sec7qrs; i++, q++)
            printf("  QRS %d: reference beat type %d\n", i, *q);

    }
    return 1;
}

