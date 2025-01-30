package main

import (
    "context"
    "github.com/go-redis/redis/v8"
     "bufio"
     "fmt"
     "io"
     "net"
     "os"
)

var ctx = context.Background()
var rdb = redis.NewClient(&redis.Options{Addr: "db.occase.de:6379", Password: "", DB: 0,})

func echo(conn net.Conn) {
	r := bufio.NewReader(conn)
	for {
		line, err := r.ReadBytes(byte('\n'))
		switch err {
		case nil:
			break
		case io.EOF:
		default:
			fmt.Println("ERROR", err)
		}

               err2 := rdb.Ping(ctx).Err()
               if err2 != nil {
                  fmt.Println("ERROR", err2)
                  panic(err2)
               }

	       conn.Write(line)
	}
}

func main() {
	l, err := net.Listen("tcp", "0.0.0.0:55555")
	if err != nil {
		fmt.Println("ERROR", err)
		os.Exit(1)
	}

	for {
		conn, err := l.Accept()
		if err != nil {
			fmt.Println("ERROR", err)
			continue
		}

		go echo(conn)
	}
}
